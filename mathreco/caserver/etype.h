#ifndef ETYPE_H_
#define ETYPE_H_

#define OP_LEAF 1
#define OP_GRP 2

#define GRPID(n) (((n) << 2) | OP_GRP)
#define LEAFID(n) (((n) << 2) | OP_LEAF)

#define OP_PARENS GRPID(1)
#define OP_SUPSCR GRPID(2)
#define OP_SUBSCR GRPID(3)
#define OP_FUNC GRPID(4)
#define OP_FRAC GRPID(5)
#define OP_MUL GRPID(6)
#define OP_ADD GRPID(7)
#define OP_SUB GRPID(8)
#define OP_EQ GRPID(9)
#define OP_NEQ GRPID(10)
#define OP_LT GRPID(11)
#define OP_LEQ GRPID(12)
#define OP_GT GRPID(13)
#define OP_GEQ GRPID(14)
#define OP_ROOT GRPID(15)
#define OP_NEG GRPID(16)
#define OP_MATRIX GRPID(17)
#define OP_MATRIXROW GRPID(18)
#define OP_INTEGRAL GRPID(19)
#define OP_SUM GRPID(20)

#define OP_NUM LEAFID(1)
#define OP_NAME LEAFID(2)
#define OP_EMPTY LEAFID(3)
#define OP_INFIN LEAFID(4)

#define ISGRP(op) ((op) & OP_GRP)
#define ISLEAF(op) ((op) & OP_LEAF)

#endif
