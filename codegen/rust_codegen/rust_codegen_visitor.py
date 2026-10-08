from codegen.skeleton.skeleton_tree import SkeletonTree
from codegen.skeleton.tree_node import SkeletonTreeNode, NoBranchDispatcherNode, DispatcherNodeBase
from codegen.skeleton.tree_node import EquationAllZeroDispatcherNode, SolvedVariableDispatcherNode
from codegen.skeleton.tree_node import ExplicitSolutionNode, PolynomialSolutionNode
from codegen.skeleton.tree_node import General6DoFNumericalReduceSolutionNode
from codegen.skeleton.tree_node import LinearSinCosType_2_SolverNode
from codegen.skeleton.tree_node import SkeletonNodeVisitor
from codegen.rust_codegen.scope_manager import RustScopeManager
from codegen.rust_codegen.rust_code_printer import YaikRustCodePrinter
from solver.solved_variable import SolutionMethod, VariableSolutionClassKey
from solver.solved_variable_impl import General6DoFNumericalReduceSolutionImpl
from typing import List, Optional, Dict, Set
import sympy as sp


class NodeIkSolutionGeneratorRust(SkeletonNodeVisitor):
    """
    Generate the body of the raw ik function in Rust. Mirrors NodeIkSolutionGeneratorCpp:
    every node processor reads the queue indices assigned to its input, solves its variable
    and forwards the new solution indices to its children. The C++ lambdas become labeled
    blocks, so ``return`` from a processor is ``break 'label``.

    The generated code expects these locals in the enclosing function:
    solution_queue: &mut SolutionQueue<..>, node_index_workspace: &mut NodeIndexWorkspace<..>,
    ik_output: &mut Vec<[f64; ROBOT_NQ]>
    """

    def __init__(self, tree: SkeletonTree, use_safe_operator: bool = False):
        super().__init__()
        self._use_safe_operator = use_safe_operator
        self._printer = YaikRustCodePrinter(use_safe_operator=use_safe_operator)
        self._tree = tree
        self._n_nodes = tree.num_nodes

        # From symbol to variable index
        self._symbol_to_index: Dict[sp.Symbol, int] = dict()
        self._index_to_symbol: Dict[int, sp.Symbol] = dict()
        self._all_variables_in_tree = tree.all_solved_symbol()

        # Build the map
        all_solved_symbol_name_list = sorted([symbol_i.name for symbol_i in self._all_variables_in_tree])
        for i in range(len(all_solved_symbol_name_list)):
            self._index_to_symbol[i] = sp.Symbol(all_solved_symbol_name_list[i])
            self._symbol_to_index[sp.Symbol(all_solved_symbol_name_list[i])] = i

        self._code_lines: List[str] = list()
        self._scope_manager = RustScopeManager(self._code_lines)

    # Code line utilities
    def add_code_line(self, code_line: str):
        self._scope_manager.append_code_line(code_line)

    def add_empty_line(self):
        self.add_code_line('')

    def enter_scope(self, scope_name: str, header: Optional[str] = None):
        self._scope_manager.enter_scope(scope_name, header)

    def exit_scope(self, scope_name: Optional[str]):
        self._scope_manager.exit_scope(scope_name)

    def chain_scope(self, scope_name: str, header: str = 'else'):
        self._scope_manager.chain_scope(scope_name, header)

    @property
    def indent_level(self) -> int:
        return self._scope_manager.indent_level

    @property
    def generated_code_lines(self):
        return self._code_lines

    @property
    def symbol_to_index(self):
        return self._symbol_to_index

    def rust_code(self, expr: sp.Expr) -> str:
        return self._printer.doprint(expr)

    def rust_name(self, symbol: sp.Symbol) -> str:
        return self._printer.doprint(symbol)

    @staticmethod
    def make_raw_solution_code() -> str:
        return 'make_intermediate_solution::<INTERMEDIATE_SOLUTION_SIZE>()'

    @staticmethod
    def add_input_index_code(node_idx: int, solution_idx_option_code: str) -> str:
        return 'node_index_workspace.add_input_index_to({node_idx}, {idx});'.format(
            node_idx=node_idx, idx=solution_idx_option_code)

    def append_init_lines(self):
        self.add_empty_line()
        self.add_code_line('solution_queue.reset();')
        self.add_code_line('node_index_workspace.reset(N_TREE_NODES);')
        assert self.indent_level == 0

    def append_finalize_lines(self):
        tree = self._tree
        robot = tree.robot
        solved_symbols = [elem.symbol for elem in tree.solved_unknowns]
        free_symbols = self.free_symbols()

        assert self.indent_level == 0
        self.add_empty_line()
        self.add_code_line('// Collect the output')
        self.enter_scope('for_all_queue', 'for i in 0..solution_queue.size()')
        self.add_code_line('if !solution_queue.solutions_validity[i] {')
        self.add_code_line('    continue;')
        self.add_code_line('}')
        self.add_code_line('let raw_ik_out_i = solution_queue.get_solution(i);')
        self.add_code_line('let mut new_ik_i = [0.0_f64; ROBOT_NQ];')
        for i in range(len(robot.unknowns)):
            unknown_i = robot.unknowns[i]
            if unknown_i.symbol in solved_symbols:
                self.add_code_line('let value_at_{idx} = raw_ik_out_i[{out_idx}];  // {symbol}'.format(
                    idx=i, out_idx=self._symbol_to_index[unknown_i.symbol], symbol=unknown_i.symbol.name))
            else:
                # This should be a parameter
                assert unknown_i.symbol in free_symbols
                self.add_code_line('let value_at_{idx} = {symbol};  // {symbol}'.format(
                    idx=i, symbol=self.rust_name(unknown_i.symbol)))
            self.add_code_line('new_ik_i[{idx}] = value_at_{idx};'.format(idx=i))
        self.add_code_line('ik_output.push(new_ik_i);')
        self.exit_scope('for_all_queue')

    # Shared pieces
    @staticmethod
    def processor_label(node: SkeletonTreeNode) -> str:
        return "'node_{idx}_processor".format(idx=node.flatten_idx)

    def _generate_processor_function_define(self, node: SkeletonTreeNode, input_idx: int):
        label = self.processor_label(node)
        self.add_code_line('// ' + self.node_processor_name(node))
        self.enter_scope(label, label + ':')
        self.add_code_line('let this_node_input_index = *node_index_workspace.get_input_indices_for_node({});'.
                           format(input_idx))
        self.add_code_line('let this_input_valid = node_index_workspace.is_input_indices_valid_for_node({});'.
                           format(input_idx))
        self.add_code_line('if !this_input_valid {')
        self.add_code_line('    break {};'.format(label))
        self.add_code_line('}')

    def _exit_processor(self, node: SkeletonTreeNode):
        self.exit_scope(self.processor_label(node))

    def _enter_for_each_input(self, solution_name: str = 'this_solution'):
        self.enter_scope('for_i', 'for &node_input_i_idx_in_queue in this_node_input_index.as_slice()')
        self.add_code_line('if !solution_queue.is_solution_valid(node_input_i_idx_in_queue) {')
        self.add_code_line('    continue;')
        self.add_code_line('}')
        self.add_code_line('let {sol_name} = *solution_queue.get_solution(node_input_i_idx_in_queue);'.format(
            sol_name=solution_name))

    def _add_temp_variables(self, replacements):
        self.add_code_line('// Temp variable for efficiency')
        for symbol_k, expr_k in replacements:
            self.add_code_line('let {var}: f64 = {code};'.format(var=self.rust_name(symbol_k),
                                                                code=self.rust_code(expr_k)))
        self.add_code_line('// End of temp variables')

    # Dispatchers
    def visit_no_branch_dispatcher(self, dispatcher: DispatcherNodeBase):
        assert dispatcher.non_degenerate_child is not None
        self.add_empty_line()
        self.add_code_line('// Code for non-branch dispatcher node {idx}'.format(idx=dispatcher.flatten_idx))
        self.add_code_line('// Actually, there is no code')

    def visit_equation_all_zero_dispatcher(self, dispatcher: EquationAllZeroDispatcherNode):
        input_idx = dispatcher.flatten_idx
        assert dispatcher.non_degenerate_child is not None
        self.add_empty_line()
        self.add_code_line('// Code for equation all-zero dispatcher node {idx}'.format(idx=input_idx))
        if dispatcher.is_root:
            self._visit_equation_all_zero_dispatcher_root(dispatcher)
        else:
            self._visit_equation_all_zero_dispatcher_non_root(dispatcher)
        self.add_code_line('// Finish code for equation all-zero dispatcher node {idx}'.format(idx=input_idx))

    def _visit_equation_all_zero_dispatcher_root(self, dispatcher: EquationAllZeroDispatcherNode):
        assert self.indent_level == 0
        checker_code = self.rust_code(dispatcher.degenerate_checker)
        non_degenerate_child_idx = dispatcher.non_degenerate_child.flatten_idx
        degenerate_child_idx = None if dispatcher.degenerate_child is None else dispatcher.degenerate_child.flatten_idx

        # The root has no input solution, nothing to extract
        self.enter_scope('checker_processor')
        self.add_code_line('let checked_result: bool = {};'.format(checker_code))
        self.enter_scope('if_non_degenerate', 'if !checked_result')
        self.add_code_line('// To non-degenerate node')
        self.add_code_line('node_index_workspace.node_input_validity_vector[{}] = true;'.format(
            non_degenerate_child_idx))
        self.exit_scope('if_non_degenerate')
        if degenerate_child_idx is not None:
            self.chain_scope('else_degenerate')
            self.add_code_line('node_index_workspace.node_input_validity_vector[{}] = true;'.format(
                degenerate_child_idx))
            self.exit_scope('else_degenerate')
        self.exit_scope('checker_processor')

    def _visit_equation_all_zero_dispatcher_non_root(self, dispatcher: EquationAllZeroDispatcherNode):
        assert self.indent_level == 0
        checker_i = dispatcher.degenerate_checker
        checker_code = self.rust_code(checker_i)
        non_degenerate_child_idx = dispatcher.non_degenerate_child.flatten_idx
        degenerate_child_idx = None if dispatcher.degenerate_child is None else dispatcher.degenerate_child.flatten_idx

        self._generate_processor_function_define(dispatcher, dispatcher.flatten_idx)
        solution_name = 'this_solution'
        requested_vars = self.get_requested_solved_variables(None, checker_i, dispatcher)
        var_extract_lines = self.extract_solved_variable_codes(requested_vars, solution_name)

        self._enter_for_each_input(solution_name)
        for this_line in var_extract_lines:
            self.add_code_line(this_line)
        self.add_code_line('let checked_result: bool = {};'.format(checker_code))
        self.enter_scope('if_non_degenerate', 'if !checked_result')
        self.add_code_line('// To non-degenerate node')
        self.add_code_line(self.add_input_index_code(non_degenerate_child_idx, 'Some(node_input_i_idx_in_queue)'))
        self.exit_scope('if_non_degenerate')
        if degenerate_child_idx is not None:
            self.chain_scope('else_degenerate')
            self.add_code_line(self.add_input_index_code(degenerate_child_idx, 'Some(node_input_i_idx_in_queue)'))
            self.exit_scope('else_degenerate')
        self.exit_scope('for_i')
        self._exit_processor(dispatcher)

    def visit_solved_variable_dispatcher(self, dispatcher: SolvedVariableDispatcherNode):
        input_idx = dispatcher.flatten_idx
        assert dispatcher.non_degenerate_child is not None
        self.add_empty_line()
        self.add_code_line('// Code for solved_variable dispatcher node {idx}'.format(idx=input_idx))
        self._generate_processor_function_define(dispatcher, input_idx)
        if dispatcher.is_root:
            raise NotImplementedError('Solved variable dispatcher can not be the root')
        self._visit_solved_variable_dispatcher_non_root(dispatcher)
        self._exit_processor(dispatcher)
        assert self.indent_level == 0
        self.add_code_line('// Finish code for solved_variable dispatcher node {idx}'.format(idx=input_idx))

    def _visit_solved_variable_dispatcher_non_root(self, dispatcher: SolvedVariableDispatcherNode):
        assert self.indent_level == 1
        n_child = len(dispatcher.branch_conditions)
        non_degenerate_child_idx = dispatcher.non_degenerate_child.flatten_idx

        # Gather all the variables used by any branch condition
        solution_name = 'this_solution'
        requested_vars = self.gather_request_solved_variables(list(dispatcher.branch_conditions), dispatcher)
        extract_lines = self.extract_solved_variable_codes(requested_vars, solution_name)

        self._enter_for_each_input(solution_name)
        self.add_code_line('let mut taken_by_degenerate = false;')
        for this_line in extract_lines:
            self.add_code_line(this_line)

        for i in range(n_child):
            child_i = dispatcher.degenerate_branch_child(i)
            condition_code = self.rust_code(dispatcher.branch_conditions[i])
            self.add_empty_line()
            self.add_code_line('let degenerate_valid_{idx}: bool = {code};'.format(idx=i, code=condition_code))
            self.enter_scope('if_degenerate', 'if degenerate_valid_{idx}'.format(idx=i))
            self.add_code_line('taken_by_degenerate = true;')
            if child_i is not None:
                self.add_code_line(self.add_input_index_code(child_i.flatten_idx, 'Some(node_input_i_idx_in_queue)'))
            self.exit_scope('if_degenerate')

        self.add_empty_line()
        self.enter_scope('if_non_degenerate', 'if !taken_by_degenerate')
        self.add_code_line(self.add_input_index_code(non_degenerate_child_idx, 'Some(node_input_i_idx_in_queue)'))
        self.exit_scope('if_non_degenerate')
        self.exit_scope('for_i')

    # Explicit solution
    def visit_explicit_solution_node(self, solution_node: ExplicitSolutionNode):
        parent_node = solution_node.parent
        assert parent_node is not None
        input_idx = parent_node.flatten_idx if isinstance(parent_node, NoBranchDispatcherNode) \
            else solution_node.flatten_idx

        self.add_empty_line()
        self.add_code_line('// Code for explicit solution node {idx}, solved variable is {var}'.
                           format(idx=solution_node.flatten_idx, var=solution_node.solution.solved_variable.name))
        self._generate_processor_function_define(solution_node, input_idx)
        if solution_node.is_parent_root_dispatcher():
            self._visit_explicit_solution_node_root(solution_node)
        else:
            self._visit_explicit_solution_node_non_root(solution_node)
        self._exit_processor(solution_node)
        self.add_code_line('// Finish code for explicit solution node {idx}'.format(idx=input_idx))
        assert self.indent_level == 0

    def _visit_explicit_solution_node_root(self, solution_node: ExplicitSolutionNode):
        assert self.indent_level == 1
        output_idx = None if solution_node.children[0] is None else solution_node.children[0].flatten_idx
        solution = solution_node.solution
        assert solution.is_explicit_solution
        solved_var_idx = self._symbol_to_index[solution.solved_variable]

        self.add_empty_line()
        self.add_code_line('// The explicit solution of root node')
        n_solution = len(solution.explicit_solutions)
        replacements, reduced_solutions = sp.cse(list(solution.explicit_solutions))
        for i in range(n_solution):
            checker_code = self.rust_code(solution.argument_validity_checkers[i])
            self.add_code_line('let condition_{i}: bool = {code};'.format(i=i, code=checker_code))
            self.enter_scope('if_condition', 'if condition_{i}'.format(i=i))
            self._add_temp_variables(replacements)
            self.add_empty_line()
            self.add_code_line('let mut solution_{i} = {make};'.format(i=i, make=self.make_raw_solution_code()))
            self.add_code_line('solution_{i}[{var_idx}] = {code};'.format(
                i=i, var_idx=solved_var_idx, code=self.rust_code(reduced_solutions[i])))
            self.add_code_line('let appended_idx = solution_queue.push_solution(solution_{i});'.format(i=i))
            if output_idx is not None:
                self.add_code_line(self.add_input_index_code(output_idx, 'appended_idx'))
            self.exit_scope('if_condition')
            self.add_empty_line()

    def _visit_explicit_solution_node_non_root(self, solution_node: ExplicitSolutionNode):
        assert self.indent_level == 1
        output_idx = None if solution_node.children[0] is None else solution_node.children[0].flatten_idx
        solution = solution_node.solution
        solved_var_idx = self._symbol_to_index[solution.solved_variable]

        self.add_empty_line()
        self.add_code_line('// The solution of non-root node {idx}'.format(idx=solution_node.flatten_idx))
        n_solution = len(solution.explicit_solutions)
        solution_list: List[sp.Expr] = list(solution.explicit_solutions)
        checker_list: List[sp.Expr] = list(solution.argument_validity_checkers)
        solution_replacements, reduced_solutions = sp.cse(solution_list)

        solution_name = 'this_solution'
        requested_vars = self.gather_request_solved_variables(solution_list + checker_list, solution_node)
        var_extract_lines = self.extract_solved_variable_codes(requested_vars, solution_name)

        self._enter_for_each_input(solution_name)
        for this_line in var_extract_lines:
            self.add_code_line(this_line)

        for i in range(n_solution):
            checker_code = self.rust_code(checker_list[i])
            solution_code = self.rust_code(reduced_solutions[i])
            self.add_empty_line()
            self.add_code_line('let condition_{i}: bool = {code};'.format(i=i, code=checker_code))
            self.enter_scope('if_condition', 'if condition_{i}'.format(i=i))

            # The temp variables may contain operations only valid when the condition holds
            self._add_temp_variables(solution_replacements)

            if i != n_solution - 1:
                self.add_code_line('let mut solution_{i} = {sol_name};'.format(i=i, sol_name=solution_name))
                self.add_code_line('solution_{i}[{var_idx}] = {code};'.format(
                    i=i, var_idx=solved_var_idx, code=solution_code))
                self.add_code_line('let appended_idx = solution_queue.push_solution(solution_{i});'.format(i=i))
                if output_idx is not None:
                    self.add_code_line(self.add_input_index_code(output_idx, 'appended_idx'))
                self.exit_scope('if_condition')
            else:
                # The last solution is written into the input solution in place
                self.add_code_line('let tmp_sol_value: f64 = {code};'.format(code=solution_code))
                self.add_code_line('solution_queue.get_solution_mut(node_input_i_idx_in_queue)[{var_idx}] '
                                   '= tmp_sol_value;'.format(var_idx=solved_var_idx))
                if output_idx is not None:
                    self.add_code_line(self.add_input_index_code(output_idx, 'Some(node_input_i_idx_in_queue)'))
                self.exit_scope('if_condition')
                self.chain_scope('else')
                self.add_code_line('solution_queue.solutions_validity[node_input_i_idx_in_queue] = false;')
                self.exit_scope('else')
        self.exit_scope('for_i')

    # Polynomial solution
    def visit_polynomial_solution_node(self, solution_node: PolynomialSolutionNode):
        parent_node = solution_node.parent
        assert parent_node is not None
        input_idx = parent_node.flatten_idx if isinstance(parent_node, NoBranchDispatcherNode) \
            else solution_node.flatten_idx

        self.add_empty_line()
        self.add_code_line('// Code for polynomial solution node {idx}, solved variable is {var}'.
                           format(idx=input_idx, var=solution_node.solution.solved_variable.name))
        self._generate_processor_function_define(solution_node, input_idx)
        if not solution_node.is_parent_root_dispatcher():
            raise NotImplementedError('We should use polynomial only if nothing-else can be found')
        self._visit_polynomial_solution_node_root(solution_node)
        self._exit_processor(solution_node)
        assert self.indent_level == 0
        self.add_code_line('// Finish code for polynomial solution node {idx}'.format(idx=input_idx))

    def _visit_polynomial_solution_node_root(self, solution_node: PolynomialSolutionNode):
        assert self.indent_level == 1
        output_idx = None if solution_node.children[0] is None else solution_node.children[0].flatten_idx
        solution = solution_node.solution
        assert solution.is_polynomial
        solved_var_idx = self._symbol_to_index[solution.solved_variable]

        self.add_empty_line()
        self.add_code_line('// The polynomial solution of root node')
        poly_dict_input = solution.polynomial_to_solve
        poly_order = -1
        for order in poly_dict_input:
            numerator, denominator = poly_dict_input[order]
            self.add_code_line('let poly_coefficient_{order}_num: f64 = {code};'.format(
                order=order, code=self.rust_code(numerator)))
            self.add_code_line('let poly_coefficient_{order}_denom: f64 = {code};'.format(
                order=order, code=self.rust_code(denominator)))
            self.add_code_line('let poly_coefficient_{order}: f64 = poly_coefficient_{order}_num / '
                               'poly_coefficient_{order}_denom;'.format(order=order))
            poly_order = max(poly_order, order)

        # Coefficients from the highest order. Orders absent from the dict are zero.
        assert poly_order >= 1
        self.add_code_line('let mut p_coefficients = [0.0_f64; {n}];'.format(n=poly_order + 1))
        for order in poly_dict_input:
            self.add_code_line('p_coefficients[{idx}] = poly_coefficient_{order};'.format(
                idx=poly_order - order, order=order))

        self.add_empty_line()
        self.add_code_line('// Invoke the solver. Note that p_coefficients[0] is the highest order')
        self.add_code_line('let poly_roots = compute_polynomial_real_roots::<{order}>(&p_coefficients);'.format(
            order=poly_order))

        self.add_empty_line()
        self.add_code_line('// Result collection')
        self.enter_scope('for_root', 'for this_root_record in poly_roots.iter()')
        self.add_code_line('if !this_root_record.is_valid {')
        self.add_code_line('    continue;')
        self.add_code_line('}')
        self.add_code_line('let this_root = this_root_record.value;')

        solution_method = solution.solution_method
        if solution_method not in (SolutionMethod.PolynomialSin.name, SolutionMethod.PolynomialCos.name):
            raise NotImplementedError('Polynomial solver codegen with non sin/cos intermediate is not implemented yet')
        self.add_code_line('if this_root.abs() > 1.0 {')
        self.add_code_line('    continue;')
        self.add_code_line('}')
        if solution_method == SolutionMethod.PolynomialSin.name:
            self.add_code_line('let first_angle = this_root.asin();')
            self.add_code_line('let second_angle = PI - this_root.asin();')
        else:
            self.add_code_line('let first_angle = this_root.acos();')
            self.add_code_line('let second_angle = -this_root.acos();')
        for sol_idx, angle in ((0, 'first_angle'), (1, 'second_angle')):
            self.add_code_line('let mut solution_{i} = {make};'.format(i=sol_idx, make=self.make_raw_solution_code()))
            self.add_code_line('solution_{i}[{var_idx}] = {angle};'.format(
                i=sol_idx, var_idx=solved_var_idx, angle=angle))
        self.add_code_line('let appended_idx_0 = solution_queue.push_solution(solution_0);')
        self.add_code_line('let appended_idx_1 = solution_queue.push_solution(solution_1);')
        if output_idx is not None:
            self.add_code_line(self.add_input_index_code(output_idx, 'appended_idx_0'))
            self.add_code_line(self.add_input_index_code(output_idx, 'appended_idx_1'))
        self.exit_scope('for_root')

    # Linear type 2 solution
    def visit_linear_type2_solution_node(self, solution_node: LinearSinCosType_2_SolverNode):
        parent_node = solution_node.parent
        assert parent_node is not None
        input_idx = parent_node.flatten_idx if isinstance(parent_node, NoBranchDispatcherNode) \
            else solution_node.flatten_idx

        self.add_empty_line()
        self.add_code_line('// Code for linear type2 solution node {idx}, solved variable is {var}'.
                           format(idx=solution_node.flatten_idx, var=solution_node.solution.solved_variable.name))
        self._generate_processor_function_define(solution_node, input_idx)
        if solution_node.is_parent_root_dispatcher():
            raise NotImplementedError('Linear type2 solution is not supported at the root')
        self._visit_linear_type2_solution_non_root(solution_node)
        self._exit_processor(solution_node)
        self.add_code_line('// Finish code for linear type2 solution node {idx}'.format(idx=input_idx))
        assert self.indent_level == 0

    def _visit_linear_type2_solution_non_root(self, solution_node: LinearSinCosType_2_SolverNode):
        assert self.indent_level == 1
        output_idx = None if solution_node.children[0] is None else solution_node.children[0].flatten_idx
        solved_var_idx = self._symbol_to_index[solution_node.solution.solved_variable]
        A_matrix: sp.Matrix = solution_node.A_matrix
        rows_to_try = solution_node.rows_to_try
        assert A_matrix.cols == 4

        self.add_empty_line()
        self.add_code_line('// The solution of non-root node {idx}'.format(idx=solution_node.flatten_idx))
        A_matrix_flat_expr: List[sp.Expr] = [A_matrix[r, c] for r in range(A_matrix.rows) for c in range(A_matrix.cols)]
        solution_replacements, reduced_A_matrix_flat = sp.cse(A_matrix_flat_expr)

        solution_name = 'this_solution'
        requested_vars = self.gather_request_solved_variables(A_matrix_flat_expr, solution_node)
        var_extract_lines = self.extract_solved_variable_codes(requested_vars, solution_name)

        self._enter_for_each_input(solution_name)
        for this_line in var_extract_lines:
            self.add_code_line(this_line)
        self._add_temp_variables(solution_replacements)

        self.add_code_line('let mut A_matrix = [[0.0_f64; 4]; {rows}];'.format(rows=A_matrix.rows))
        offset = 0
        for r in range(A_matrix.rows):
            for c in range(A_matrix.cols):
                rc_expr: sp.Expr = reduced_A_matrix_flat[offset]
                if rc_expr != sp.S.Zero:
                    self.add_code_line('A_matrix[{r}][{c}] = {code};'.format(r=r, c=c, code=self.rust_code(rc_expr)))
                offset += 1

        self.add_code_line('// Code for solving')
        for rows_j in rows_to_try:
            r0, r1, r2 = rows_j
            header = ('if let Some((var_solution_0, var_solution_1)) = '
                      'linear_solver::try_solve_linear_type2_specific_rows::<{rows}>(&A_matrix, {r0}, {r1}, {r2})'.
                      format(rows=A_matrix.rows, r0=r0, r1=r1, r2=r2))
            self.enter_scope('if_solved', header)
            # The first candidate is a new solution
            self.add_code_line('let mut solution_0 = {sol_name};'.format(sol_name=solution_name))
            self.add_code_line('solution_0[{var_idx}] = var_solution_0;'.format(var_idx=solved_var_idx))
            self.add_code_line('let appended_idx = solution_queue.push_solution(solution_0);')
            if output_idx is not None:
                self.add_code_line(self.add_input_index_code(output_idx, 'appended_idx'))
            # The second candidate is written in place
            self.add_code_line('solution_queue.get_solution_mut(node_input_i_idx_in_queue)[{var_idx}] = '
                               'var_solution_1;'.format(var_idx=solved_var_idx))
            if output_idx is not None:
                self.add_code_line(self.add_input_index_code(output_idx, 'Some(node_input_i_idx_in_queue)'))
            self.add_code_line('continue;')
            self.exit_scope('if_solved')
        self.exit_scope('for_i')

    # General 6-dof solution
    def visit_general_6dof_numerical_reduce_node(self, solution_node: General6DoFNumericalReduceSolutionNode):
        parent_node = solution_node.parent
        assert parent_node is not None
        input_idx = parent_node.flatten_idx if isinstance(parent_node, NoBranchDispatcherNode) \
            else solution_node.flatten_idx

        self.add_empty_line()
        self.add_code_line('// Code for general_6dof solution node {idx}, solved variable is {var}'.
                           format(idx=solution_node.flatten_idx, var=solution_node.solution.solved_variable.name))
        self._generate_processor_function_define(solution_node, input_idx)
        if not solution_node.is_parent_root_dispatcher():
            raise NotImplementedError('General 6-DoF solution is only supported at the root')
        if solution_node.has_semi_symbolic_reduce:
            self._visit_general_6dof_semi_symbolic_reduce_node_root(solution_node)
        else:
            self._visit_general_6dof_numerical_reduce_node_root(solution_node)
        self._exit_processor(solution_node)
        self.add_code_line('// Finish code for general_6dof solution node {idx}'.format(idx=input_idx))
        assert self.indent_level == 0

    def _generate_dmatrix_code(self, matrix_name: str, n_rows: int, n_cols: int, entries: List[sp.Expr]):
        self.add_code_line('let mut {name} = DMatrix::<f64>::zeros({r}, {c});'.format(
            name=matrix_name, r=n_rows, c=n_cols))
        for r in range(n_rows):
            for c in range(n_cols):
                rc_expr = entries[r * n_cols + c]
                if rc_expr != sp.S.Zero:
                    self.add_code_line('{name}[({r}, {c})] = {code};'.format(
                        name=matrix_name, r=r, c=c, code=self.rust_code(rc_expr)))

    def _generate_dmatrix_code_from_sp(self, sp_matrix: sp.Matrix, matrix_name: str):
        n_rows, n_cols = sp_matrix.shape[0], sp_matrix.shape[1]
        entries = [sp_matrix[r, c] for r in range(n_rows) for c in range(n_cols)]
        self._generate_dmatrix_code(matrix_name, n_rows, n_cols, entries)

    def _append_general_6dof_solutions(self, solution_node: General6DoFNumericalReduceSolutionNode):
        output_idx = None if solution_node.children[0] is None else solution_node.children[0].flatten_idx
        solved_var_idx = self._symbol_to_index[solution_node.solved_variable]
        self.add_empty_line()
        self.enter_scope('for_loop', 'for &solution_value in solution_buffer[..n_solutions].iter()')
        self.add_code_line('let mut solution_i = {make};'.format(make=self.make_raw_solution_code()))
        self.add_code_line('solution_i[{var_idx}] = solution_value;'.format(var_idx=solved_var_idx))
        self.add_code_line('let appended_idx = solution_queue.push_solution(solution_i);')
        if output_idx is not None:
            self.add_code_line(self.add_input_index_code(output_idx, 'appended_idx'))
        self.exit_scope('for_loop')

    def _visit_general_6dof_semi_symbolic_reduce_node_root(self, solution_node: General6DoFNumericalReduceSolutionNode):
        assert self.indent_level == 1
        assert solution_node.has_semi_symbolic_reduce
        solution_impl = solution_node.general_6dof_solution.internal_solution
        assert solution_impl.is_general_6dof
        assert solution_impl.class_key() == VariableSolutionClassKey.General6DoFNumericalReduce.name
        general_6dof_solution: General6DoFNumericalReduceSolutionImpl = solution_impl
        semi_reduce = general_6dof_solution.semi_symbolic_reduce_record
        label = self.processor_label(solution_node)

        self.add_empty_line()
        self.add_code_line('// The general 6-dof solution of root node')
        R_l = semi_reduce.R_l
        self._generate_dmatrix_code_from_sp(R_l, 'R_l')
        self.add_code_line('let R_l_mat_inv = match R_l.try_inverse() {')
        self.add_code_line('    Some(inv) if general_6dof_internal::all_finite(&inv) => inv,')
        self.add_code_line('    _ => break {},'.format(label))
        self.add_code_line('};')

        self.add_empty_line()
        R_l_inv_symbols = semi_reduce.R_l_inv_as_symbols
        for r in range(R_l.shape[0]):
            for c in range(R_l.shape[1]):
                self.add_code_line('let {name}: f64 = R_l_mat_inv[({r}, {c})];'.format(
                    name=self.rust_name(R_l_inv_symbols[r, c]), r=r, c=c))

        A, B, C = semi_reduce.A, semi_reduce.B, semi_reduce.C
        matrix_expressions: List[sp.Expr] = list()
        for sp_matrix in (A, B, C):
            matrix_expressions.extend(sp_matrix[r, c] for r in range(sp_matrix.shape[0]) for c in range(sp_matrix.shape[1]))
        replacements, reduced_expressions = sp.cse(matrix_expressions)

        self.add_empty_line()
        self._add_temp_variables(replacements)
        counter = 0
        for sp_matrix, name in ((A, 'A'), (B, 'B'), (C, 'C')):
            n_rows, n_cols = sp_matrix.shape[0], sp_matrix.shape[1]
            self.add_empty_line()
            self._generate_dmatrix_code(name, n_rows, n_cols, reduced_expressions[counter:counter + n_rows * n_cols])
            counter += n_rows * n_cols

        self.add_empty_line()
        self.add_code_line('// Invoke the solver')
        self.add_code_line('let mut solution_buffer = [0.0_f64; 16];')
        self.add_code_line('let n_solutions = general_6dof_internal::compute_solution_from_tanhalf_lme'
                           '(&A, &B, &C, &mut solution_buffer);')
        self._append_general_6dof_solutions(solution_node)

    def _visit_general_6dof_numerical_reduce_node_root(self, solution_node: General6DoFNumericalReduceSolutionNode):
        assert self.indent_level == 1
        solution_impl = solution_node.general_6dof_solution.internal_solution
        assert solution_impl.is_general_6dof
        assert solution_impl.class_key() == VariableSolutionClassKey.General6DoFNumericalReduce.name
        general_6dof_solution: General6DoFNumericalReduceSolutionImpl = solution_impl

        self.add_empty_line()
        self.add_code_line('// The general 6-dof solution of root node')
        A_sin, A_cos, C_const = general_6dof_solution.lhs_matrices()
        N_rhs = general_6dof_solution.rhs_matrix()
        for sp_matrix, name in ((A_sin, 'A_sin'), (A_cos, 'A_cos'), (C_const, 'C_const'), (N_rhs, 'N_rhs')):
            self._generate_dmatrix_code_from_sp(sp_matrix, name)
            self.add_empty_line()

        n_total_lines = A_sin.shape[0]
        lines2reduce: List[int] = list(general_6dof_solution.select_lines)
        remaining_lines: List[int] = [i for i in range(n_total_lines) if i not in lines2reduce]

        self.add_code_line('// The lines for reduction')
        self.add_code_line('let lines_to_reduce: [usize; {n}] = [{v}];'.format(
            n=len(lines2reduce), v=', '.join(str(elem) for elem in lines2reduce)))
        self.add_code_line('let remaining_lines: [usize; {n}] = [{v}];'.format(
            n=len(remaining_lines), v=', '.join(str(elem) for elem in remaining_lines)))

        self.add_empty_line()
        self.add_code_line('// Invoke the solver')
        self.add_code_line('let mut solution_buffer = [0.0_f64; 16];')
        self.add_code_line('let n_solutions = general_6dof_numerical_reduce_solve(&A_sin, &A_cos, &C_const, &N_rhs, '
                           '&lines_to_reduce, &remaining_lines, &mut solution_buffer);')
        self._append_general_6dof_solutions(solution_node)

    # Naming and variable utilities
    @staticmethod
    def node_processor_name(node: SkeletonTreeNode) -> str:
        type_name = type(node).__name__
        if node.is_dispatcher_node:
            return '{type_name}_node_{idx}_processor'.format(idx=node.flatten_idx, type_name=type_name)
        solved_var_name: str = node.solution.solved_variable.name
        return '{type_name}_node_{idx}_solve_{var_name}_processor'.format(
            idx=node.flatten_idx, type_name=type_name, var_name=solved_var_name)

    def gather_request_solved_variables(self, expressions_in: List[sp.Expr], node: SkeletonTreeNode) -> List[sp.Symbol]:
        requested_symbol_names_set: Set[str] = set()
        for expr_i in expressions_in:
            for symbol_j in expr_i.free_symbols:
                if symbol_j in self._all_variables_in_tree:
                    requested_symbol_names_set.add(symbol_j.name)
        requested_symbol_names = sorted(requested_symbol_names_set)
        assert self._verify_requested_variables_solved(requested_symbol_names, node), \
            'Node {} requests variables not solved by its parents: {}'.format(node.flatten_idx, requested_symbol_names)
        return [sp.Symbol(symbol_name) for symbol_name in requested_symbol_names]

    def _verify_requested_variables_solved(self, requested_symbol_names: List[str], node: SkeletonTreeNode) -> bool:
        solutions_from_root = self._tree.variable_solutions_from_root(node, include_current_node=False)
        solved_names = set(solution_j.solved_variable.name for solution_j in solutions_from_root)
        return all(name in solved_names for name in requested_symbol_names)

    def extract_solved_variable_codes(self, requested_variables: List[sp.Symbol],
                                      solution_name: str = 'this_solution') -> List[str]:
        return ['let {var_name}: f64 = {sol_name}[{idx}];'.format(
            var_name=self.rust_name(var_i), sol_name=solution_name, idx=self._symbol_to_index[var_i])
            for var_i in requested_variables]

    def get_requested_solved_variables(self, solution_i: Optional[sp.Expr], checker_i: Optional[sp.Expr],
                                       node: SkeletonTreeNode) -> List[sp.Symbol]:
        expressions = [expr for expr in (solution_i, checker_i) if expr is not None]
        return self.gather_request_solved_variables(expressions, node)

    def free_symbols(self, include_unknown_as_parameter: bool = True) -> List[sp.Symbol]:
        return self._tree.get_free_symbols(include_unknown_as_parameter)
