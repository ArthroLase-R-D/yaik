from codegen.skeleton.codegen_base import indent
from typing import List, Optional


class RustScopeManager(object):
    """
    Track brace-delimited scopes and indentation for the generated Rust code.
    A scope may carry a header, e.g. ``if condition`` or a block label ``'node_3_processor:``.
    """

    def __init__(self, code_lines_workspace: List[str]):
        self._code_lines = code_lines_workspace
        self._scope_name_stack: List[str] = list()

    def enter_scope(self, scope_name: str, header: Optional[str] = None):
        if header is None:
            self._code_lines.append(self.current_indent() + '{')
        else:
            self._code_lines.append(self.current_indent() + header + ' {')
        self._scope_name_stack.append(scope_name)

    def exit_scope(self, scope_name: Optional[str], trailing: str = ''):
        if scope_name is not None:
            assert self._scope_name_stack[-1] == scope_name, \
                'Expect to exit scope {}, but the top scope is {}'.format(scope_name, self._scope_name_stack[-1])
        self._scope_name_stack.pop()
        self._code_lines.append(self.current_indent() + '}' + trailing)

    def chain_scope(self, scope_name: str, header: str = 'else'):
        """
        Continue the scope that was just closed, e.g. turn the closing ``}`` of an if into ``} else {``.
        """
        assert self._code_lines[-1] == self.current_indent() + '}'
        self._code_lines[-1] = self.current_indent() + '} ' + header + ' {'
        self._scope_name_stack.append(scope_name)

    def current_indent(self) -> str:
        return indent(self.indent_level)

    def append_code_line(self, code_line: str):
        if len(code_line) == 0:
            self._code_lines.append('')
        else:
            self._code_lines.append(self.current_indent() + code_line)

    @property
    def indent_level(self):
        return len(self._scope_name_stack)
