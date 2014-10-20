# This Python script parses the following C++ header files 
# where constants are defined that need to be referenced in SageWrapper:
#    include/grammar-values.h
#    caserver/cmdcode.h

import re # regular expressions

# This method removes C++ comments
# @Author: Markus Jarderot
def remove_comments(text):
    def replacer(match):
        s = match.group(0)
        if s.startswith('/'):
            return ""
        else:
            return s
    pattern = re.compile(
        r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',
        re.DOTALL | re.MULTILINE
    )
    return re.sub(pattern, replacer, text)

expr_list = []
cmd_list = []

file = open("../include/grammar-values.h", "r").read()
lines = remove_comments(file).splitlines()

suffix = '_EXPR'
r = 'SemanticId +\w+%s *= *-? *[0-9]+ *;' % suffix

for line in lines:
    m = re.search(r, line)
    if m:   
        s = m.group(0)
        m = re.search('\w+%s' % suffix, s)
        name = m.group(0)
        m = re.search('-? *[0-9]+', s)
        value = m.group(0)
        expr_list.append(name)
        exec("%s = %s" % (name, value));

file = open("cmdcode.h", "r").read()
lines = remove_comments(file).splitlines()

prefix = 'CASCMD_'
r = '#define +%s\w+ +-? *[0-9]+' % prefix

for line in lines:
    m = re.search(r, line)
    if m:   
        s = m.group(0)
        m = re.search('%s\w+' % prefix, s)
        name = m.group(0)
        m = re.search('-? *[0-9]+', s)
        value = m.group(0)
        cmd_list.append(name)
        exec("%s = %s" % (name, value));

# Returns a convenient string identifying the expression type
def GetExprString(code):
    
    # An exception
    if code == TERMINAL_EXPR:
        return "TERM"
    
    for var in expr_list:
        if globals()[var] == code:
            return var[:len(var)-len(suffix)] # remove the suffix
        
    return repr(code) # UNRECOGNIZED

