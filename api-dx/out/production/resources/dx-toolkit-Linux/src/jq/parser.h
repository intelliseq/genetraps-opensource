/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_YY_PARSER_H_INCLUDED
# define YY_YY_PARSER_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif
/* "%code requires" blocks.  */
#line 9 "parser.y" /* yacc.c:1909  */

#include "locfile.h"
struct lexer_param;

#define YYLTYPE location
#define YYLLOC_DEFAULT(Loc, Rhs, N)             \
  do {                                          \
    if (N) {                                    \
      (Loc).start = YYRHSLOC(Rhs, 1).start;     \
      (Loc).end = YYRHSLOC(Rhs, N).end;         \
    } else {                                    \
      (Loc).start = YYRHSLOC(Rhs, 0).end;       \
      (Loc).end = YYRHSLOC(Rhs, 0).end;         \
    }                                           \
  } while (0)
 

#line 62 "parser.h" /* yacc.c:1909  */

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    INVALID_CHARACTER = 258,
    IDENT = 259,
    FIELD = 260,
    LITERAL = 261,
    FORMAT = 262,
    SETMOD = 263,
    EQ = 264,
    NEQ = 265,
    DEFINEDOR = 266,
    AS = 267,
    DEF = 268,
    IF = 269,
    THEN = 270,
    ELSE = 271,
    ELSE_IF = 272,
    REDUCE = 273,
    END = 274,
    AND = 275,
    OR = 276,
    SETPIPE = 277,
    SETPLUS = 278,
    SETMINUS = 279,
    SETMULT = 280,
    SETDIV = 281,
    SETDEFINEDOR = 282,
    LESSEQ = 283,
    GREATEREQ = 284,
    QQSTRING_START = 285,
    QQSTRING_TEXT = 286,
    QQSTRING_INTERP_START = 287,
    QQSTRING_INTERP_END = 288,
    QQSTRING_END = 289
  };
#endif
/* Tokens.  */
#define INVALID_CHARACTER 258
#define IDENT 259
#define FIELD 260
#define LITERAL 261
#define FORMAT 262
#define SETMOD 263
#define EQ 264
#define NEQ 265
#define DEFINEDOR 266
#define AS 267
#define DEF 268
#define IF 269
#define THEN 270
#define ELSE 271
#define ELSE_IF 272
#define REDUCE 273
#define END 274
#define AND 275
#define OR 276
#define SETPIPE 277
#define SETPLUS 278
#define SETMINUS 279
#define SETMULT 280
#define SETDIV 281
#define SETDEFINEDOR 282
#define LESSEQ 283
#define GREATEREQ 284
#define QQSTRING_START 285
#define QQSTRING_TEXT 286
#define QQSTRING_INTERP_START 287
#define QQSTRING_INTERP_END 288
#define QQSTRING_END 289

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 29 "parser.y" /* yacc.c:1909  */

  jv literal;
  block blk;

#line 147 "parser.h" /* yacc.c:1909  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif

/* Location type.  */
#if ! defined YYLTYPE && ! defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE YYLTYPE;
struct YYLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
};
# define YYLTYPE_IS_DECLARED 1
# define YYLTYPE_IS_TRIVIAL 1
#endif



int yyparse (block* answer, int* errors, struct locfile* locations, struct lexer_param* lexer_param_ptr);

#endif /* !YY_YY_PARSER_H_INCLUDED  */
