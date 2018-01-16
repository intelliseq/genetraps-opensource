/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison implementation for Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 1 "parser.y" /* yacc.c:339  */

#include <stdio.h>
#include <string.h>
#include "compile.h"
#include "jv_alloc.h"
#define YYMALLOC jv_mem_alloc
#define YYFREE jv_mem_free

#line 75 "parser.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 1
#endif

/* In a future release of Bison, this section will be replaced
   by #include "y.tab.h".  */
#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif
/* "%code requires" blocks.  */
#line 9 "parser.y" /* yacc.c:355  */

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
 

#line 123 "parser.c" /* yacc.c:355  */

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
#line 29 "parser.y" /* yacc.c:355  */

  jv literal;
  block blk;

#line 208 "parser.c" /* yacc.c:355  */
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

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */
#line 95 "parser.y" /* yacc.c:358  */

#include "lexer.h"
struct lexer_param {
  yyscan_t lexer;
};
#define FAIL(loc, msg)                                             \
  do {                                                             \
    location l = loc;                                              \
    yyerror(&l, answer, errors, locations, lexer_param_ptr, msg);  \
    /*YYERROR*/;                                                   \
  } while (0)

void yyerror(YYLTYPE* loc, block* answer, int* errors, 
             struct locfile* locations, struct lexer_param* lexer_param_ptr, const char *s){
  (*errors)++;
  locfile_locate(locations, *loc, "error: %s", s);
}

int yylex(YYSTYPE* yylval, YYLTYPE* yylloc, block* answer, int* errors, 
          struct locfile* locations, struct lexer_param* lexer_param_ptr) {
  yyscan_t lexer = lexer_param_ptr->lexer;
  int tok = jq_yylex(yylval, yylloc, lexer);
  if ((tok == LITERAL || tok == QQSTRING_TEXT) && !jv_is_valid(yylval->literal)) {
    jv msg = jv_invalid_get_msg(jv_copy(yylval->literal));
    if (jv_get_kind(msg) == JV_KIND_STRING) {
      FAIL(*yylloc, jv_string_value(msg));
    } else {
      FAIL(*yylloc, "Invalid literal");
    }
    jv_free(msg);
    jv_free(yylval->literal);
    yylval->literal = jv_null();
  }
  return tok;
}

static block gen_dictpair(block k, block v) {
  return BLOCK(gen_subexp(k), gen_subexp(v), gen_op_simple(INSERT));
}

static block gen_index(block obj, block key) {
  return BLOCK(gen_subexp(key), obj, gen_op_simple(INDEX));
}

static block gen_slice_index(block obj, block start, block end) {
  block key = BLOCK(gen_subexp(gen_const(jv_object())),
                    gen_subexp(gen_const(jv_string("start"))),
                    gen_subexp(start),
                    gen_op_simple(INSERT),
                    gen_subexp(gen_const(jv_string("end"))),
                    gen_subexp(end),
                    gen_op_simple(INSERT));
  return BLOCK(key, obj, gen_op_simple(INDEX));
}

static block gen_binop(block a, block b, int op) {
  const char* funcname = 0;
  switch (op) {
  case '+': funcname = "_plus"; break;
  case '-': funcname = "_minus"; break;
  case '*': funcname = "_multiply"; break;
  case '/': funcname = "_divide"; break;
  case '%': funcname = "_mod"; break;
  case EQ: funcname = "_equal"; break;
  case NEQ: funcname = "_notequal"; break;
  case '<': funcname = "_less"; break;
  case '>': funcname = "_greater"; break;
  case LESSEQ: funcname = "_lesseq"; break;
  case GREATEREQ: funcname = "_greatereq"; break;
  }
  assert(funcname);

  return gen_call(funcname, BLOCK(gen_lambda(a), gen_lambda(b)));
}

static block gen_format(block a, jv fmt) {
  return BLOCK(a, gen_call("format", BLOCK(gen_lambda(gen_const(fmt)))));
}

static block gen_definedor_assign(block object, block val) {
  block tmp = gen_op_var_fresh(STOREV, "tmp");
  return BLOCK(gen_op_simple(DUP),
               val, tmp,
               gen_call("_modify", BLOCK(gen_lambda(object),
                                         gen_lambda(gen_definedor(gen_noop(), 
                                                                  gen_op_bound(LOADV, tmp))))));
}
 
static block gen_update(block object, block val, int optype) {
  block tmp = gen_op_var_fresh(STOREV, "tmp");
  return BLOCK(gen_op_simple(DUP),
               val,
               tmp,
               gen_call("_modify", BLOCK(gen_lambda(object), 
                                         gen_lambda(gen_binop(gen_noop(),
                                                              gen_op_bound(LOADV, tmp),
                                                              optype)))));
}


#line 336 "parser.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL \
             && defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
  YYLTYPE yyls_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE) + sizeof (YYLTYPE)) \
      + 2 * YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  45
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1294

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  55
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  14
/* YYNRULES -- Number of rules.  */
#define YYNRULES  87
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  193

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   289

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,    46,    45,     2,     2,
      47,    48,    43,    41,    37,    42,    50,    44,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    49,    35,
      39,    38,    40,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    51,     2,    52,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    53,    36,    54,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   198,   198,   201,   206,   209,   214,   218,   223,   228,
     231,   236,   240,   244,   248,   252,   256,   260,   264,   268,
     272,   276,   280,   284,   288,   292,   296,   300,   304,   308,
     312,   316,   320,   324,   328,   332,   336,   341,   346,   354,
     367,   367,   371,   371,   378,   381,   384,   390,   393,   398,
     402,   408,   411,   414,   417,   420,   423,   427,   433,   436,
     439,   442,   445,   448,   451,   454,   457,   460,   463,   466,
     469,   473,   477,   482,   487,   488,   489,   490,   493,   496,
     497,   498,   501,   504,   507,   511,   515,   518
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "INVALID_CHARACTER", "IDENT", "FIELD",
  "LITERAL", "FORMAT", "\"%=\"", "\"==\"", "\"!=\"", "\"//\"", "\"as\"",
  "\"def\"", "\"if\"", "\"then\"", "\"else\"", "\"elif\"", "\"reduce\"",
  "\"end\"", "\"and\"", "\"or\"", "\"|=\"", "\"+=\"", "\"-=\"", "\"*=\"",
  "\"/=\"", "\"//=\"", "\"<=\"", "\">=\"", "QQSTRING_START",
  "QQSTRING_TEXT", "QQSTRING_INTERP_START", "QQSTRING_INTERP_END",
  "QQSTRING_END", "';'", "'|'", "','", "'='", "'<'", "'>'", "'+'", "'-'",
  "'*'", "'/'", "'%'", "'$'", "'('", "')'", "':'", "'.'", "'['", "']'",
  "'{'", "'}'", "$accept", "TopLevel", "FuncDefs", "Exp", "FuncDef",
  "String", "@1", "@2", "QQString", "ElseBody", "ExpD", "Term", "MkDict",
  "MkDictPair", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,    59,   124,    44,    61,    60,
      62,    43,    45,    42,    47,    37,    36,    40,    41,    58,
      46,    91,    93,   123,   125
};
# endif

#define YYPACT_NINF -89

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-89)))

#define YYTABLE_NINF -79

#define yytable_value_is_error(Yytable_value) \
  (!!((Yytable_value) == (-79)))

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     348,   -20,   -89,   -89,   -15,    30,   348,   398,   -89,   348,
      33,   252,   185,   237,    82,    39,   -89,  1138,   348,   -89,
      26,   348,   -89,   -21,   834,   348,    31,   -89,   -26,   -89,
     -18,   632,   -89,    41,   -15,   -89,   -12,   -89,   477,   -25,
      -5,   267,    29,    34,     8,   -89,   348,   348,   348,   348,
     348,   348,   348,   348,   348,   348,   348,   348,   348,   348,
     348,   348,   348,   348,   348,   348,   348,   348,   348,   348,
     -89,  1138,   -89,    47,    -6,   100,   591,   -89,    81,   348,
     318,    48,    -9,   -89,   -89,   -89,   -89,   -89,    91,   -89,
     398,    49,   673,   398,   -89,    91,  1202,   378,   378,  1176,
    1249,  1228,  1202,  1202,  1202,  1202,  1202,  1202,   378,   378,
    1138,  1176,  1202,   378,   378,   -26,   -26,   -89,   -89,   -89,
      92,   -89,    50,   348,   -89,   432,   348,   -89,     1,   -32,
     872,   -89,   796,    95,   -89,   348,   -89,    63,   -89,    67,
      66,    59,    60,    67,   -89,    74,   -89,   515,   333,   -89,
     714,   -89,   107,    70,   -89,   348,   348,   -89,    68,   910,
     398,   398,   398,   348,   -89,   -89,   553,   -89,    72,   348,
     948,   986,   348,   -89,    67,    67,    67,  1138,   -89,    73,
    1024,   -89,   348,  1062,   348,   -89,   796,   348,  1100,   -89,
     755,   -89,   -89
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       4,    71,    53,    63,    65,     0,     0,     0,    40,     0,
       0,     0,     0,     0,     0,     0,     3,     2,     4,    64,
      36,     0,    42,     0,     0,     0,     0,    44,    21,    70,
       0,     0,    56,     0,     0,    55,     0,    68,     0,     0,
      85,     0,    84,     0,    79,     1,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       5,     6,    52,     0,     0,     0,     0,    44,     0,     0,
       0,     0,     0,    74,    66,    57,    75,    67,     0,    77,
       0,     0,     0,     0,    69,     0,    29,    30,    31,    14,
      13,    12,    16,    20,    23,    25,    28,    15,    34,    35,
      17,    18,    11,    32,    33,    19,    22,    24,    26,    27,
       0,    54,     0,     0,    59,     0,     0,    72,     0,     0,
       0,    10,     0,     0,    45,     0,    41,     0,    81,    82,
      50,     0,     0,    83,    80,     0,    76,     0,     0,    58,
       0,    43,     0,     0,    37,     0,     0,     9,     0,     0,
       0,     0,     0,     0,    62,    61,     0,    73,     0,     0,
       0,     0,     0,    46,    49,    87,    86,     7,    60,     0,
       0,    48,     0,     0,     0,    38,     0,     0,     0,    47,
       0,    39,     8
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -89,   -89,   106,     0,     2,    -4,   -89,   -89,    51,   -61,
     -88,    -3,   -81,   -89
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    15,    16,    71,    25,    19,    27,    77,    82,   157,
     139,    20,    43,    44
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      17,    34,    18,   152,    26,   143,    24,   138,    35,    28,
      42,    31,    88,    38,   144,    22,   153,    67,    68,    69,
      18,    76,   134,   135,     8,   136,    78,    21,    79,    89,
      83,    72,   134,   135,    23,   151,    72,    29,    73,    45,
      86,    92,    85,    81,    90,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,   108,   109,
     110,   111,   112,   113,   114,   115,   116,   117,   118,   119,
     121,    72,   174,   175,   176,   125,    74,    75,    93,   130,
     132,    74,    75,    39,    42,   129,    40,   140,    94,    34,
     140,    42,   137,   120,   133,    40,   145,   141,    34,   158,
      88,   122,   146,   160,     1,     2,     3,     4,   161,   162,
     163,   168,     8,     5,     6,   172,    74,    75,     7,   169,
     179,     8,   184,   147,    70,   189,   150,     0,   128,    41,
       8,     0,     0,     0,     0,   159,   -78,     0,    41,     0,
       0,     0,     9,     0,     0,   -78,    10,    11,   166,   123,
      12,    13,   124,    14,     0,   170,   171,   140,   140,   140,
       0,     0,     0,   177,     0,     0,     0,     0,     0,   180,
       0,     0,   183,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   186,     0,   188,   -51,    32,   190,     0,    33,
     -51,     0,    34,   -51,   -51,   -51,   -51,   -51,     0,     0,
     -51,   -51,   -51,     0,   -51,   -51,   -51,   -51,   -51,   -51,
     -51,   -51,   -51,   -51,   -51,     8,     0,     0,   -51,     0,
     -51,   -51,   -51,   -51,   -51,   -51,   -51,   -51,   -51,   -51,
     -51,     0,     0,   -51,   -51,   -51,   -51,   -51,    36,   -51,
       0,     1,     2,     3,     4,     0,     0,     0,     0,     0,
       5,     6,     0,    30,     0,     7,     1,     2,     3,     4,
       0,     0,     0,     0,     0,     5,     6,     8,    91,     0,
       7,     1,     2,     3,     4,     0,     0,     0,     0,     9,
       5,     6,     8,    10,    11,     7,     0,    12,    13,    37,
      14,     0,     0,     0,     9,     0,     0,     8,    10,    11,
       0,     0,    12,    13,     0,    14,     0,     0,     0,     9,
       0,     0,     0,    10,    11,     0,     0,    12,    13,   131,
      14,     0,     1,     2,     3,     4,     0,     0,     0,     0,
       0,     5,     6,     0,     0,     0,     7,     1,     2,     3,
       4,     0,     0,     0,     0,     0,     5,     6,     8,     0,
       0,     7,     1,     2,     3,     4,     0,     0,     0,     0,
       9,     5,     6,     8,    10,    11,     7,     0,    12,    13,
       0,    14,     0,     0,     0,     9,     0,     0,     8,    10,
      11,     0,     0,    12,    13,   165,    14,   -79,   -79,     0,
       9,     0,     0,     0,    10,    11,     0,     0,    12,    13,
       0,    14,     1,     2,     3,     4,   -79,   -79,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   -79,   -79,    65,
      66,    67,    68,    69,     0,     0,     0,     0,     8,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      46,    47,    48,    49,    10,    11,     0,     0,    12,    13,
       0,    14,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,     0,     0,     0,     0,     0,     0,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    69,     0,     0,
       0,   148,     0,     0,   149,    46,    47,    48,    49,     0,
       0,     0,     0,     0,     0,     0,     0,    50,    51,    52,
      53,    54,    55,    56,    57,    58,    59,     0,     0,     0,
       0,     0,     0,    60,    61,    62,    63,    64,    65,    66,
      67,    68,    69,    46,    47,    48,    49,     0,     0,    87,
       0,     0,     0,     0,     0,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,     0,     0,     0,     0,     0,
       0,    60,    61,    62,    63,    64,    65,    66,    67,    68,
      69,    46,    47,    48,    49,     0,     0,   164,     0,     0,
       0,     0,     0,    50,    51,    52,    53,    54,    55,    56,
      57,    58,    59,     0,     0,     0,     0,     0,     0,    60,
      61,    62,    63,    64,    65,    66,    67,    68,    69,    46,
      47,    48,    49,     0,     0,   178,     0,     0,     0,     0,
       0,    50,    51,    52,    53,    54,    55,    56,    57,    58,
      59,     0,     0,     0,     0,     0,   126,    60,    61,    62,
      63,    64,    65,    66,    67,    68,    69,     0,     0,   127,
      46,    47,    48,    49,     0,     0,     0,     0,     0,     0,
       0,     0,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,     0,     0,     0,     0,     0,     0,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    69,     0,     0,
      84,    46,    47,    48,    49,     0,     0,     0,     0,     0,
       0,     0,     0,    50,    51,    52,    53,    54,    55,    56,
      57,    58,    59,     0,     0,     0,     0,     0,     0,    60,
      61,    62,    63,    64,    65,    66,    67,    68,    69,     0,
       0,   142,    46,    47,    48,    49,     0,     0,     0,     0,
       0,     0,     0,     0,    50,    51,    52,    53,    54,    55,
      56,    57,    58,    59,     0,     0,     0,     0,     0,     0,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    69,
       0,     0,   167,    46,    47,    48,    49,     0,     0,     0,
       0,     0,     0,     0,     0,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,     0,     0,     0,     0,     0,
       0,    60,    61,    62,    63,    64,    65,    66,    67,    68,
      69,     0,     0,   192,    46,    47,    48,    49,     0,     0,
       0,     0,   155,   156,     0,     0,    50,    51,    52,    53,
      54,    55,    56,    57,    58,    59,     0,     0,     0,     0,
       0,     0,    60,    61,    62,    63,    64,    65,    66,    67,
      68,    69,    46,    47,    48,    49,     0,     0,     0,    80,
       0,     0,     0,     0,    50,    51,    52,    53,    54,    55,
      56,    57,    58,    59,     0,     0,     0,     0,     0,     0,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    69,
      46,    47,    48,    49,     0,     0,     0,     0,     0,     0,
       0,     0,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,     0,     0,     0,     0,     0,   154,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    69,    46,    47,
      48,    49,     0,     0,     0,     0,     0,     0,     0,     0,
      50,    51,    52,    53,    54,    55,    56,    57,    58,    59,
       0,     0,     0,   173,     0,     0,    60,    61,    62,    63,
      64,    65,    66,    67,    68,    69,    46,    47,    48,    49,
       0,     0,     0,     0,     0,     0,     0,   181,    50,    51,
      52,    53,    54,    55,    56,    57,    58,    59,     0,     0,
       0,     0,     0,     0,    60,    61,    62,    63,    64,    65,
      66,    67,    68,    69,    46,    47,    48,    49,     0,     0,
       0,   182,     0,     0,     0,     0,    50,    51,    52,    53,
      54,    55,    56,    57,    58,    59,     0,     0,     0,     0,
       0,     0,    60,    61,    62,    63,    64,    65,    66,    67,
      68,    69,    46,    47,    48,    49,     0,     0,     0,     0,
       0,     0,     0,     0,    50,    51,    52,    53,    54,    55,
      56,    57,    58,    59,     0,     0,     0,     0,     0,   185,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    69,
      46,    47,    48,    49,     0,     0,     0,     0,     0,     0,
       0,     0,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,     0,     0,     0,     0,     0,   187,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    69,    46,    47,
      48,    49,     0,     0,     0,     0,     0,     0,     0,     0,
      50,    51,    52,    53,    54,    55,    56,    57,    58,    59,
       0,     0,     0,     0,     0,   191,    60,    61,    62,    63,
      64,    65,    66,    67,    68,    69,    46,    47,    48,    49,
       0,     0,     0,     0,     0,     0,     0,     0,    50,    51,
      52,    53,    54,    55,    56,    57,    58,    59,     0,     0,
       0,     0,     0,     0,    60,    61,    62,    63,    64,    65,
      66,    67,    68,    69,    46,    47,    48,    49,     0,     0,
       0,     0,     0,     0,     0,     0,    50,    51,    52,    53,
      54,    55,    56,    57,    58,    59,     0,     0,     0,     0,
     -79,    47,    48,     0,    62,    63,    64,    65,    66,    67,
      68,    69,    50,    51,   -79,   -79,   -79,   -79,   -79,   -79,
      58,    59,     0,     0,     0,     0,     0,    47,    48,     0,
     -79,    63,    64,    65,    66,    67,    68,    69,    50,     0,
       0,     0,     0,     0,     0,     0,    58,    59,    47,    48,
       0,     0,     0,     0,     0,     0,     0,    63,    64,    65,
      66,    67,    68,    69,     0,     0,     0,    58,    59,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    63,    64,
      65,    66,    67,    68,    69
};

static const yytype_int16 yycheck[] =
{
       0,     7,     0,    35,     7,    93,     6,    88,    12,     9,
      14,    11,    37,    13,    95,    30,    48,    43,    44,    45,
      18,    21,    31,    32,    30,    34,    47,    47,    49,    54,
      48,     5,    31,    32,     4,    34,     5,     4,    12,     0,
      52,    41,     1,    12,    49,    37,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    56,    57,    58,    59,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    69,
      74,     5,   160,   161,   162,    75,    50,    51,    49,    79,
      80,    50,    51,     1,    88,     4,     4,    90,    54,     7,
      93,    95,     1,    46,    46,     4,     4,    48,     7,     4,
      37,     1,    52,    36,     4,     5,     6,     7,    49,    49,
      36,     4,    30,    13,    14,    47,    50,    51,    18,    49,
      48,    30,    49,   123,    18,   186,   126,    -1,    77,    47,
      30,    -1,    -1,    -1,    -1,   135,    54,    -1,    47,    -1,
      -1,    -1,    42,    -1,    -1,    54,    46,    47,   148,    49,
      50,    51,    52,    53,    -1,   155,   156,   160,   161,   162,
      -1,    -1,    -1,   163,    -1,    -1,    -1,    -1,    -1,   169,
      -1,    -1,   172,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   182,    -1,   184,     0,     1,   187,    -1,     4,
       5,    -1,     7,     8,     9,    10,    11,    12,    -1,    -1,
      15,    16,    17,    -1,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    -1,    -1,    33,    -1,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    -1,    -1,    48,    49,    50,    51,    52,     1,    54,
      -1,     4,     5,     6,     7,    -1,    -1,    -1,    -1,    -1,
      13,    14,    -1,     1,    -1,    18,     4,     5,     6,     7,
      -1,    -1,    -1,    -1,    -1,    13,    14,    30,     1,    -1,
      18,     4,     5,     6,     7,    -1,    -1,    -1,    -1,    42,
      13,    14,    30,    46,    47,    18,    -1,    50,    51,    52,
      53,    -1,    -1,    -1,    42,    -1,    -1,    30,    46,    47,
      -1,    -1,    50,    51,    -1,    53,    -1,    -1,    -1,    42,
      -1,    -1,    -1,    46,    47,    -1,    -1,    50,    51,     1,
      53,    -1,     4,     5,     6,     7,    -1,    -1,    -1,    -1,
      -1,    13,    14,    -1,    -1,    -1,    18,     4,     5,     6,
       7,    -1,    -1,    -1,    -1,    -1,    13,    14,    30,    -1,
      -1,    18,     4,     5,     6,     7,    -1,    -1,    -1,    -1,
      42,    13,    14,    30,    46,    47,    18,    -1,    50,    51,
      -1,    53,    -1,    -1,    -1,    42,    -1,    -1,    30,    46,
      47,    -1,    -1,    50,    51,    52,    53,     9,    10,    -1,
      42,    -1,    -1,    -1,    46,    47,    -1,    -1,    50,    51,
      -1,    53,     4,     5,     6,     7,    28,    29,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    39,    40,    41,
      42,    43,    44,    45,    -1,    -1,    -1,    -1,    30,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
       8,     9,    10,    11,    46,    47,    -1,    -1,    50,    51,
      -1,    53,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    -1,    -1,    -1,    -1,    -1,    -1,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    -1,    -1,
      -1,    49,    -1,    -1,    52,     8,     9,    10,    11,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    -1,    -1,    -1,
      -1,    -1,    -1,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,     8,     9,    10,    11,    -1,    -1,    52,
      -1,    -1,    -1,    -1,    -1,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    -1,    -1,    -1,    -1,    -1,
      -1,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,     8,     9,    10,    11,    -1,    -1,    52,    -1,    -1,
      -1,    -1,    -1,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    -1,    -1,    -1,    -1,    -1,    -1,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,     8,
       9,    10,    11,    -1,    -1,    52,    -1,    -1,    -1,    -1,
      -1,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    -1,    -1,    -1,    -1,    -1,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    -1,    -1,    48,
       8,     9,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    -1,    -1,    -1,    -1,    -1,    -1,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    -1,    -1,
      48,     8,     9,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    -1,    -1,    -1,    -1,    -1,    -1,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    -1,
      -1,    48,     8,     9,    10,    11,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      -1,    -1,    48,     8,     9,    10,    11,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    -1,    -1,    -1,    -1,    -1,
      -1,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    -1,    -1,    48,     8,     9,    10,    11,    -1,    -1,
      -1,    -1,    16,    17,    -1,    -1,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    -1,    -1,    -1,    -1,
      -1,    -1,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,     8,     9,    10,    11,    -1,    -1,    -1,    15,
      -1,    -1,    -1,    -1,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    -1,    -1,    -1,    -1,    -1,    -1,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
       8,     9,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    -1,    -1,    -1,    -1,    -1,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,     8,     9,
      10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      -1,    -1,    -1,    33,    -1,    -1,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,     8,     9,    10,    11,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    -1,    -1,
      -1,    -1,    -1,    -1,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,     8,     9,    10,    11,    -1,    -1,
      -1,    15,    -1,    -1,    -1,    -1,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    -1,    -1,    -1,    -1,
      -1,    -1,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,     8,     9,    10,    11,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    -1,    -1,    -1,    -1,    -1,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
       8,     9,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    -1,    -1,    -1,    -1,    -1,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,     8,     9,
      10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      -1,    -1,    -1,    -1,    -1,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,     8,     9,    10,    11,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    -1,    -1,
      -1,    -1,    -1,    -1,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,     8,     9,    10,    11,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    -1,    -1,    -1,    -1,
       8,     9,    10,    -1,    38,    39,    40,    41,    42,    43,
      44,    45,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    -1,    -1,    -1,    -1,    -1,     9,    10,    -1,
      38,    39,    40,    41,    42,    43,    44,    45,    20,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    28,    29,     9,    10,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    39,    40,    41,
      42,    43,    44,    45,    -1,    -1,    -1,    28,    29,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    39,    40,
      41,    42,    43,    44,    45
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     4,     5,     6,     7,    13,    14,    18,    30,    42,
      46,    47,    50,    51,    53,    56,    57,    58,    59,    60,
      66,    47,    30,     4,    58,    59,    66,    61,    58,     4,
       1,    58,     1,     4,     7,    60,     1,    52,    58,     1,
       4,    47,    60,    67,    68,     0,     8,     9,    10,    11,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      57,    58,     5,    12,    50,    51,    58,    62,    47,    49,
      15,    12,    63,    48,    48,     1,    52,    52,    37,    54,
      49,     1,    58,    49,    54,    37,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      46,    60,     1,    49,    52,    58,    35,    48,    63,     4,
      58,     1,    58,    46,    31,    32,    34,     1,    67,    65,
      66,    48,    48,    65,    67,     4,    52,    58,    49,    52,
      58,    34,    35,    48,    35,    16,    17,    64,     4,    58,
      36,    49,    49,    36,    52,    52,    58,    48,     4,    49,
      58,    58,    47,    33,    65,    65,    65,    58,    52,    48,
      58,    19,    15,    58,    49,    35,    58,    35,    58,    64,
      58,    35,    48
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    55,    56,    56,    57,    57,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    58,    58,    58,    58,    58,    59,    59,    59,
      61,    60,    62,    60,    63,    63,    63,    64,    64,    65,
      65,    66,    66,    66,    66,    66,    66,    66,    66,    66,
      66,    66,    66,    66,    66,    66,    66,    66,    66,    66,
      66,    66,    66,    66,    66,    66,    66,    66,    67,    67,
      67,    67,    68,    68,    68,    68,    68,    68
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     0,     2,     2,     6,    10,     5,
       4,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     2,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     1,     5,     8,    10,
       0,     4,     0,     5,     0,     2,     4,     5,     3,     3,
       1,     1,     2,     1,     3,     2,     2,     3,     4,     3,
       6,     5,     5,     1,     1,     1,     3,     3,     2,     3,
       2,     1,     4,     6,     3,     3,     4,     3,     0,     1,
       3,     3,     3,     3,     1,     1,     5,     5
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (&yylloc, answer, errors, locations, lexer_param_ptr, YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)                                \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;        \
          (Current).first_column = YYRHSLOC (Rhs, 1).first_column;      \
          (Current).last_line    = YYRHSLOC (Rhs, N).last_line;         \
          (Current).last_column  = YYRHSLOC (Rhs, N).last_column;       \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).first_line   = (Current).last_line   =              \
            YYRHSLOC (Rhs, 0).last_line;                                \
          (Current).first_column = (Current).last_column =              \
            YYRHSLOC (Rhs, 0).last_column;                              \
        }                                                               \
    while (0)
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K])


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL

/* Print *YYLOCP on YYO.  Private, do not rely on its existence. */

YY_ATTRIBUTE_UNUSED
static unsigned
yy_location_print_ (FILE *yyo, YYLTYPE const * const yylocp)
{
  unsigned res = 0;
  int end_col = 0 != yylocp->last_column ? yylocp->last_column - 1 : 0;
  if (0 <= yylocp->first_line)
    {
      res += YYFPRINTF (yyo, "%d", yylocp->first_line);
      if (0 <= yylocp->first_column)
        res += YYFPRINTF (yyo, ".%d", yylocp->first_column);
    }
  if (0 <= yylocp->last_line)
    {
      if (yylocp->first_line < yylocp->last_line)
        {
          res += YYFPRINTF (yyo, "-%d", yylocp->last_line);
          if (0 <= end_col)
            res += YYFPRINTF (yyo, ".%d", end_col);
        }
      else if (0 <= end_col && yylocp->first_column < end_col)
        res += YYFPRINTF (yyo, "-%d", end_col);
    }
  return res;
 }

#  define YY_LOCATION_PRINT(File, Loc)          \
  yy_location_print_ (File, &(Loc))

# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, Location, answer, errors, locations, lexer_param_ptr); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, block* answer, int* errors, struct locfile* locations, struct lexer_param* lexer_param_ptr)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (yylocationp);
  YYUSE (answer);
  YYUSE (errors);
  YYUSE (locations);
  YYUSE (lexer_param_ptr);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, block* answer, int* errors, struct locfile* locations, struct lexer_param* lexer_param_ptr)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  YY_LOCATION_PRINT (yyoutput, *yylocationp);
  YYFPRINTF (yyoutput, ": ");
  yy_symbol_value_print (yyoutput, yytype, yyvaluep, yylocationp, answer, errors, locations, lexer_param_ptr);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, YYLTYPE *yylsp, int yyrule, block* answer, int* errors, struct locfile* locations, struct lexer_param* lexer_param_ptr)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                       , &(yylsp[(yyi + 1) - (yynrhs)])                       , answer, errors, locations, lexer_param_ptr);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, yylsp, Rule, answer, errors, locations, lexer_param_ptr); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, YYLTYPE *yylocationp, block* answer, int* errors, struct locfile* locations, struct lexer_param* lexer_param_ptr)
{
  YYUSE (yyvaluep);
  YYUSE (yylocationp);
  YYUSE (answer);
  YYUSE (errors);
  YYUSE (locations);
  YYUSE (lexer_param_ptr);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  switch (yytype)
    {
          case 4: /* IDENT  */
#line 34 "parser.y" /* yacc.c:1257  */
      { jv_free(((*yyvaluep).literal)); }
#line 1585 "parser.c" /* yacc.c:1257  */
        break;

    case 5: /* FIELD  */
#line 34 "parser.y" /* yacc.c:1257  */
      { jv_free(((*yyvaluep).literal)); }
#line 1591 "parser.c" /* yacc.c:1257  */
        break;

    case 6: /* LITERAL  */
#line 34 "parser.y" /* yacc.c:1257  */
      { jv_free(((*yyvaluep).literal)); }
#line 1597 "parser.c" /* yacc.c:1257  */
        break;

    case 7: /* FORMAT  */
#line 34 "parser.y" /* yacc.c:1257  */
      { jv_free(((*yyvaluep).literal)); }
#line 1603 "parser.c" /* yacc.c:1257  */
        break;

    case 31: /* QQSTRING_TEXT  */
#line 34 "parser.y" /* yacc.c:1257  */
      { jv_free(((*yyvaluep).literal)); }
#line 1609 "parser.c" /* yacc.c:1257  */
        break;

    case 57: /* FuncDefs  */
#line 35 "parser.y" /* yacc.c:1257  */
      { block_free(((*yyvaluep).blk)); }
#line 1615 "parser.c" /* yacc.c:1257  */
        break;

    case 58: /* Exp  */
#line 35 "parser.y" /* yacc.c:1257  */
      { block_free(((*yyvaluep).blk)); }
#line 1621 "parser.c" /* yacc.c:1257  */
        break;

    case 59: /* FuncDef  */
#line 35 "parser.y" /* yacc.c:1257  */
      { block_free(((*yyvaluep).blk)); }
#line 1627 "parser.c" /* yacc.c:1257  */
        break;

    case 60: /* String  */
#line 35 "parser.y" /* yacc.c:1257  */
      { block_free(((*yyvaluep).blk)); }
#line 1633 "parser.c" /* yacc.c:1257  */
        break;

    case 63: /* QQString  */
#line 35 "parser.y" /* yacc.c:1257  */
      { block_free(((*yyvaluep).blk)); }
#line 1639 "parser.c" /* yacc.c:1257  */
        break;

    case 64: /* ElseBody  */
#line 35 "parser.y" /* yacc.c:1257  */
      { block_free(((*yyvaluep).blk)); }
#line 1645 "parser.c" /* yacc.c:1257  */
        break;

    case 65: /* ExpD  */
#line 35 "parser.y" /* yacc.c:1257  */
      { block_free(((*yyvaluep).blk)); }
#line 1651 "parser.c" /* yacc.c:1257  */
        break;

    case 66: /* Term  */
#line 35 "parser.y" /* yacc.c:1257  */
      { block_free(((*yyvaluep).blk)); }
#line 1657 "parser.c" /* yacc.c:1257  */
        break;

    case 67: /* MkDict  */
#line 35 "parser.y" /* yacc.c:1257  */
      { block_free(((*yyvaluep).blk)); }
#line 1663 "parser.c" /* yacc.c:1257  */
        break;

    case 68: /* MkDictPair  */
#line 35 "parser.y" /* yacc.c:1257  */
      { block_free(((*yyvaluep).blk)); }
#line 1669 "parser.c" /* yacc.c:1257  */
        break;


      default:
        break;
    }
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/*----------.
| yyparse.  |
`----------*/

int
yyparse (block* answer, int* errors, struct locfile* locations, struct lexer_param* lexer_param_ptr)
{
/* The lookahead symbol.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

/* Location data for the lookahead symbol.  */
static YYLTYPE yyloc_default
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
  = { 1, 1, 1, 1 }
# endif
;
YYLTYPE yylloc = yyloc_default;

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.
       'yyls': related to locations.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    /* The location stack.  */
    YYLTYPE yylsa[YYINITDEPTH];
    YYLTYPE *yyls;
    YYLTYPE *yylsp;

    /* The locations where the error started and ended.  */
    YYLTYPE yyerror_range[3];

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
  YYLTYPE yyloc;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N), yylsp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yylsp = yyls = yylsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  yylsp[0] = yylloc;
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;
        YYLTYPE *yyls1 = yyls;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yyls1, yysize * sizeof (*yylsp),
                    &yystacksize);

        yyls = yyls1;
        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
        YYSTACK_RELOCATE (yyls_alloc, yyls);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
      yylsp = yyls + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex (&yylval, &yylloc, answer, errors, locations, lexer_param_ptr);
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END
  *++yylsp = yylloc;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

  /* Default location.  */
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 198 "parser.y" /* yacc.c:1646  */
    {
  *answer = (yyvsp[0].blk);
}
#line 1965 "parser.c" /* yacc.c:1646  */
    break;

  case 3:
#line 201 "parser.y" /* yacc.c:1646  */
    {
  *answer = (yyvsp[0].blk);
}
#line 1973 "parser.c" /* yacc.c:1646  */
    break;

  case 4:
#line 206 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_noop();
}
#line 1981 "parser.c" /* yacc.c:1646  */
    break;

  case 5:
#line 209 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = block_join((yyvsp[-1].blk), (yyvsp[0].blk));
}
#line 1989 "parser.c" /* yacc.c:1646  */
    break;

  case 6:
#line 214 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = block_bind((yyvsp[-1].blk), (yyvsp[0].blk), OP_IS_CALL_PSEUDO);
}
#line 1997 "parser.c" /* yacc.c:1646  */
    break;

  case 7:
#line 218 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_var_binding((yyvsp[-5].blk), jv_string_value((yyvsp[-2].literal)), (yyvsp[0].blk));
  jv_free((yyvsp[-2].literal));
}
#line 2006 "parser.c" /* yacc.c:1646  */
    break;

  case 8:
#line 223 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_reduce(jv_string_value((yyvsp[-5].literal)), (yyvsp[-8].blk), (yyvsp[-3].blk), (yyvsp[-1].blk));
  jv_free((yyvsp[-5].literal));
}
#line 2015 "parser.c" /* yacc.c:1646  */
    break;

  case 9:
#line 228 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_cond((yyvsp[-3].blk), (yyvsp[-1].blk), (yyvsp[0].blk));
}
#line 2023 "parser.c" /* yacc.c:1646  */
    break;

  case 10:
#line 231 "parser.y" /* yacc.c:1646  */
    {
  FAIL((yyloc), "Possibly unterminated 'if' statment");
  (yyval.blk) = (yyvsp[-2].blk);
}
#line 2032 "parser.c" /* yacc.c:1646  */
    break;

  case 11:
#line 236 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_call("_assign", BLOCK(gen_lambda((yyvsp[-2].blk)), gen_lambda((yyvsp[0].blk))));
}
#line 2040 "parser.c" /* yacc.c:1646  */
    break;

  case 12:
#line 240 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_or((yyvsp[-2].blk), (yyvsp[0].blk));
}
#line 2048 "parser.c" /* yacc.c:1646  */
    break;

  case 13:
#line 244 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_and((yyvsp[-2].blk), (yyvsp[0].blk));
}
#line 2056 "parser.c" /* yacc.c:1646  */
    break;

  case 14:
#line 248 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_definedor((yyvsp[-2].blk), (yyvsp[0].blk));
}
#line 2064 "parser.c" /* yacc.c:1646  */
    break;

  case 15:
#line 252 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_definedor_assign((yyvsp[-2].blk), (yyvsp[0].blk));
}
#line 2072 "parser.c" /* yacc.c:1646  */
    break;

  case 16:
#line 256 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_call("_modify", BLOCK(gen_lambda((yyvsp[-2].blk)), gen_lambda((yyvsp[0].blk))));
}
#line 2080 "parser.c" /* yacc.c:1646  */
    break;

  case 17:
#line 260 "parser.y" /* yacc.c:1646  */
    { 
  (yyval.blk) = block_join((yyvsp[-2].blk), (yyvsp[0].blk)); 
}
#line 2088 "parser.c" /* yacc.c:1646  */
    break;

  case 18:
#line 264 "parser.y" /* yacc.c:1646  */
    { 
  (yyval.blk) = gen_both((yyvsp[-2].blk), (yyvsp[0].blk)); 
}
#line 2096 "parser.c" /* yacc.c:1646  */
    break;

  case 19:
#line 268 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-2].blk), (yyvsp[0].blk), '+');
}
#line 2104 "parser.c" /* yacc.c:1646  */
    break;

  case 20:
#line 272 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_update((yyvsp[-2].blk), (yyvsp[0].blk), '+');
}
#line 2112 "parser.c" /* yacc.c:1646  */
    break;

  case 21:
#line 276 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = BLOCK((yyvsp[0].blk), gen_call("_negate", gen_noop()));
}
#line 2120 "parser.c" /* yacc.c:1646  */
    break;

  case 22:
#line 280 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-2].blk), (yyvsp[0].blk), '-');
}
#line 2128 "parser.c" /* yacc.c:1646  */
    break;

  case 23:
#line 284 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_update((yyvsp[-2].blk), (yyvsp[0].blk), '-');
}
#line 2136 "parser.c" /* yacc.c:1646  */
    break;

  case 24:
#line 288 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-2].blk), (yyvsp[0].blk), '*');
}
#line 2144 "parser.c" /* yacc.c:1646  */
    break;

  case 25:
#line 292 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_update((yyvsp[-2].blk), (yyvsp[0].blk), '*');
}
#line 2152 "parser.c" /* yacc.c:1646  */
    break;

  case 26:
#line 296 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-2].blk), (yyvsp[0].blk), '/');
}
#line 2160 "parser.c" /* yacc.c:1646  */
    break;

  case 27:
#line 300 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-2].blk), (yyvsp[0].blk), '%');
}
#line 2168 "parser.c" /* yacc.c:1646  */
    break;

  case 28:
#line 304 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_update((yyvsp[-2].blk), (yyvsp[0].blk), '/');
}
#line 2176 "parser.c" /* yacc.c:1646  */
    break;

  case 29:
#line 308 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_update((yyvsp[-2].blk), (yyvsp[0].blk), '%');
}
#line 2184 "parser.c" /* yacc.c:1646  */
    break;

  case 30:
#line 312 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-2].blk), (yyvsp[0].blk), EQ);
}
#line 2192 "parser.c" /* yacc.c:1646  */
    break;

  case 31:
#line 316 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-2].blk), (yyvsp[0].blk), NEQ);
}
#line 2200 "parser.c" /* yacc.c:1646  */
    break;

  case 32:
#line 320 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-2].blk), (yyvsp[0].blk), '<');
}
#line 2208 "parser.c" /* yacc.c:1646  */
    break;

  case 33:
#line 324 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-2].blk), (yyvsp[0].blk), '>');
}
#line 2216 "parser.c" /* yacc.c:1646  */
    break;

  case 34:
#line 328 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-2].blk), (yyvsp[0].blk), LESSEQ);
}
#line 2224 "parser.c" /* yacc.c:1646  */
    break;

  case 35:
#line 332 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-2].blk), (yyvsp[0].blk), GREATEREQ);
}
#line 2232 "parser.c" /* yacc.c:1646  */
    break;

  case 36:
#line 336 "parser.y" /* yacc.c:1646  */
    { 
  (yyval.blk) = (yyvsp[0].blk); 
}
#line 2240 "parser.c" /* yacc.c:1646  */
    break;

  case 37:
#line 341 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_function(jv_string_value((yyvsp[-3].literal)), gen_noop(), (yyvsp[-1].blk));
  jv_free((yyvsp[-3].literal));
}
#line 2249 "parser.c" /* yacc.c:1646  */
    break;

  case 38:
#line 346 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_function(jv_string_value((yyvsp[-6].literal)), 
                    gen_param(jv_string_value((yyvsp[-4].literal))), 
                    (yyvsp[-1].blk));
  jv_free((yyvsp[-6].literal));
  jv_free((yyvsp[-4].literal));
}
#line 2261 "parser.c" /* yacc.c:1646  */
    break;

  case 39:
#line 354 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_function(jv_string_value((yyvsp[-8].literal)), 
                    BLOCK(gen_param(jv_string_value((yyvsp[-6].literal))), 
                          gen_param(jv_string_value((yyvsp[-4].literal)))),
                    (yyvsp[-1].blk));
  jv_free((yyvsp[-8].literal));
  jv_free((yyvsp[-6].literal));
  jv_free((yyvsp[-4].literal));
}
#line 2275 "parser.c" /* yacc.c:1646  */
    break;

  case 40:
#line 367 "parser.y" /* yacc.c:1646  */
    { (yyval.literal) = jv_string("text"); }
#line 2281 "parser.c" /* yacc.c:1646  */
    break;

  case 41:
#line 367 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = (yyvsp[-1].blk);
  jv_free((yyvsp[-2].literal));
}
#line 2290 "parser.c" /* yacc.c:1646  */
    break;

  case 42:
#line 371 "parser.y" /* yacc.c:1646  */
    { (yyval.literal) = (yyvsp[-1].literal); }
#line 2296 "parser.c" /* yacc.c:1646  */
    break;

  case 43:
#line 371 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = (yyvsp[-1].blk);
  jv_free((yyvsp[-2].literal));
}
#line 2305 "parser.c" /* yacc.c:1646  */
    break;

  case 44:
#line 378 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_const(jv_string(""));
}
#line 2313 "parser.c" /* yacc.c:1646  */
    break;

  case 45:
#line 381 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-1].blk), gen_const((yyvsp[0].literal)), '+');
}
#line 2321 "parser.c" /* yacc.c:1646  */
    break;

  case 46:
#line 384 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_binop((yyvsp[-3].blk), gen_format((yyvsp[-1].blk), jv_copy((yyvsp[-4].literal))), '+');
}
#line 2329 "parser.c" /* yacc.c:1646  */
    break;

  case 47:
#line 390 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_cond((yyvsp[-3].blk), (yyvsp[-1].blk), (yyvsp[0].blk));
}
#line 2337 "parser.c" /* yacc.c:1646  */
    break;

  case 48:
#line 393 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = (yyvsp[-1].blk);
}
#line 2345 "parser.c" /* yacc.c:1646  */
    break;

  case 49:
#line 398 "parser.y" /* yacc.c:1646  */
    { 
  (yyval.blk) = block_join((yyvsp[-2].blk), (yyvsp[0].blk));
}
#line 2353 "parser.c" /* yacc.c:1646  */
    break;

  case 50:
#line 402 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = (yyvsp[0].blk);
}
#line 2361 "parser.c" /* yacc.c:1646  */
    break;

  case 51:
#line 408 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_noop(); 
}
#line 2369 "parser.c" /* yacc.c:1646  */
    break;

  case 52:
#line 411 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_index((yyvsp[-1].blk), gen_const((yyvsp[0].literal)));
}
#line 2377 "parser.c" /* yacc.c:1646  */
    break;

  case 53:
#line 414 "parser.y" /* yacc.c:1646  */
    { 
  (yyval.blk) = gen_index(gen_noop(), gen_const((yyvsp[0].literal))); 
}
#line 2385 "parser.c" /* yacc.c:1646  */
    break;

  case 54:
#line 417 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_index((yyvsp[-2].blk), (yyvsp[0].blk));
}
#line 2393 "parser.c" /* yacc.c:1646  */
    break;

  case 55:
#line 420 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_index(gen_noop(), (yyvsp[0].blk));
}
#line 2401 "parser.c" /* yacc.c:1646  */
    break;

  case 56:
#line 423 "parser.y" /* yacc.c:1646  */
    {
  FAIL((yyloc), "try .[\"field\"] instead of .field for unusually named fields");
  (yyval.blk) = gen_noop();
}
#line 2410 "parser.c" /* yacc.c:1646  */
    break;

  case 57:
#line 427 "parser.y" /* yacc.c:1646  */
    {
  jv_free((yyvsp[-1].literal));
  FAIL((yyloc), "try .[\"field\"] instead of .field for unusually named fields");
  (yyval.blk) = gen_noop();
}
#line 2420 "parser.c" /* yacc.c:1646  */
    break;

  case 58:
#line 433 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_index((yyvsp[-3].blk), (yyvsp[-1].blk)); 
}
#line 2428 "parser.c" /* yacc.c:1646  */
    break;

  case 59:
#line 436 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = block_join((yyvsp[-2].blk), gen_op_simple(EACH)); 
}
#line 2436 "parser.c" /* yacc.c:1646  */
    break;

  case 60:
#line 439 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_slice_index((yyvsp[-5].blk), (yyvsp[-3].blk), (yyvsp[-1].blk));
}
#line 2444 "parser.c" /* yacc.c:1646  */
    break;

  case 61:
#line 442 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_slice_index((yyvsp[-4].blk), (yyvsp[-2].blk), gen_const(jv_null()));
}
#line 2452 "parser.c" /* yacc.c:1646  */
    break;

  case 62:
#line 445 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_slice_index((yyvsp[-4].blk), gen_const(jv_null()), (yyvsp[-1].blk));
}
#line 2460 "parser.c" /* yacc.c:1646  */
    break;

  case 63:
#line 448 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_const((yyvsp[0].literal)); 
}
#line 2468 "parser.c" /* yacc.c:1646  */
    break;

  case 64:
#line 451 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = (yyvsp[0].blk);
}
#line 2476 "parser.c" /* yacc.c:1646  */
    break;

  case 65:
#line 454 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_format(gen_noop(), (yyvsp[0].literal));
}
#line 2484 "parser.c" /* yacc.c:1646  */
    break;

  case 66:
#line 457 "parser.y" /* yacc.c:1646  */
    { 
  (yyval.blk) = (yyvsp[-1].blk); 
}
#line 2492 "parser.c" /* yacc.c:1646  */
    break;

  case 67:
#line 460 "parser.y" /* yacc.c:1646  */
    { 
  (yyval.blk) = gen_collect((yyvsp[-1].blk)); 
}
#line 2500 "parser.c" /* yacc.c:1646  */
    break;

  case 68:
#line 463 "parser.y" /* yacc.c:1646  */
    { 
  (yyval.blk) = gen_const(jv_array()); 
}
#line 2508 "parser.c" /* yacc.c:1646  */
    break;

  case 69:
#line 466 "parser.y" /* yacc.c:1646  */
    { 
  (yyval.blk) = BLOCK(gen_subexp(gen_const(jv_object())), (yyvsp[-1].blk), gen_op_simple(POP));
}
#line 2516 "parser.c" /* yacc.c:1646  */
    break;

  case 70:
#line 469 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_location((yyloc), gen_op_unbound(LOADV, jv_string_value((yyvsp[0].literal))));
  jv_free((yyvsp[0].literal));
}
#line 2525 "parser.c" /* yacc.c:1646  */
    break;

  case 71:
#line 473 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_location((yyloc), gen_call(jv_string_value((yyvsp[0].literal)), gen_noop()));
  jv_free((yyvsp[0].literal));
}
#line 2534 "parser.c" /* yacc.c:1646  */
    break;

  case 72:
#line 477 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_call(jv_string_value((yyvsp[-3].literal)), gen_lambda((yyvsp[-1].blk)));
  (yyval.blk) = gen_location((yylsp[-3]), (yyval.blk));
  jv_free((yyvsp[-3].literal));
}
#line 2544 "parser.c" /* yacc.c:1646  */
    break;

  case 73:
#line 482 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_call(jv_string_value((yyvsp[-5].literal)), BLOCK(gen_lambda((yyvsp[-3].blk)), gen_lambda((yyvsp[-1].blk))));
  (yyval.blk) = gen_location((yylsp[-5]), (yyval.blk));
  jv_free((yyvsp[-5].literal));
}
#line 2554 "parser.c" /* yacc.c:1646  */
    break;

  case 74:
#line 487 "parser.y" /* yacc.c:1646  */
    { (yyval.blk) = gen_noop(); }
#line 2560 "parser.c" /* yacc.c:1646  */
    break;

  case 75:
#line 488 "parser.y" /* yacc.c:1646  */
    { (yyval.blk) = gen_noop(); }
#line 2566 "parser.c" /* yacc.c:1646  */
    break;

  case 76:
#line 489 "parser.y" /* yacc.c:1646  */
    { (yyval.blk) = (yyvsp[-3].blk); }
#line 2572 "parser.c" /* yacc.c:1646  */
    break;

  case 77:
#line 490 "parser.y" /* yacc.c:1646  */
    { (yyval.blk) = gen_noop(); }
#line 2578 "parser.c" /* yacc.c:1646  */
    break;

  case 78:
#line 493 "parser.y" /* yacc.c:1646  */
    { 
  (yyval.blk)=gen_noop(); 
}
#line 2586 "parser.c" /* yacc.c:1646  */
    break;

  case 79:
#line 496 "parser.y" /* yacc.c:1646  */
    { (yyval.blk) = (yyvsp[0].blk); }
#line 2592 "parser.c" /* yacc.c:1646  */
    break;

  case 80:
#line 497 "parser.y" /* yacc.c:1646  */
    { (yyval.blk)=block_join((yyvsp[-2].blk), (yyvsp[0].blk)); }
#line 2598 "parser.c" /* yacc.c:1646  */
    break;

  case 81:
#line 498 "parser.y" /* yacc.c:1646  */
    { (yyval.blk) = (yyvsp[0].blk); }
#line 2604 "parser.c" /* yacc.c:1646  */
    break;

  case 82:
#line 501 "parser.y" /* yacc.c:1646  */
    { 
  (yyval.blk) = gen_dictpair(gen_const((yyvsp[-2].literal)), (yyvsp[0].blk));
 }
#line 2612 "parser.c" /* yacc.c:1646  */
    break;

  case 83:
#line 504 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_dictpair((yyvsp[-2].blk), (yyvsp[0].blk));
  }
#line 2620 "parser.c" /* yacc.c:1646  */
    break;

  case 84:
#line 507 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_dictpair((yyvsp[0].blk), BLOCK(gen_op_simple(POP), gen_op_simple(DUP2),
                              gen_op_simple(DUP2), gen_op_simple(INDEX)));
  }
#line 2629 "parser.c" /* yacc.c:1646  */
    break;

  case 85:
#line 511 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_dictpair(gen_const(jv_copy((yyvsp[0].literal))),
                    gen_index(gen_noop(), gen_const((yyvsp[0].literal))));
  }
#line 2638 "parser.c" /* yacc.c:1646  */
    break;

  case 86:
#line 515 "parser.y" /* yacc.c:1646  */
    {
  (yyval.blk) = gen_dictpair((yyvsp[-3].blk), (yyvsp[0].blk));
  }
#line 2646 "parser.c" /* yacc.c:1646  */
    break;

  case 87:
#line 518 "parser.y" /* yacc.c:1646  */
    { (yyval.blk) = (yyvsp[0].blk); }
#line 2652 "parser.c" /* yacc.c:1646  */
    break;


#line 2656 "parser.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;
  *++yylsp = yyloc;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (&yylloc, answer, errors, locations, lexer_param_ptr, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (&yylloc, answer, errors, locations, lexer_param_ptr, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }

  yyerror_range[1] = yylloc;

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, &yylloc, answer, errors, locations, lexer_param_ptr);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  yyerror_range[1] = yylsp[1-yylen];
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;

      yyerror_range[1] = *yylsp;
      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp, yylsp, answer, errors, locations, lexer_param_ptr);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  yyerror_range[2] = yylloc;
  /* Using YYLLOC is tempting, but would change the location of
     the lookahead.  YYLOC is available though.  */
  YYLLOC_DEFAULT (yyloc, yyerror_range, 2);
  *++yylsp = yyloc;

  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (&yylloc, answer, errors, locations, lexer_param_ptr, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, &yylloc, answer, errors, locations, lexer_param_ptr);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, yylsp, answer, errors, locations, lexer_param_ptr);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 519 "parser.y" /* yacc.c:1906  */


int jq_parse(struct locfile* locations, block* answer) {
  struct lexer_param scanner;
  YY_BUFFER_STATE buf;
  jq_yylex_init_extra(0, &scanner.lexer);
  buf = jq_yy_scan_bytes(locations->data, locations->length, scanner.lexer);
  int errors = 0;
  *answer = gen_noop();
  yyparse(answer, &errors, locations, &scanner);
  jq_yy_delete_buffer(buf, scanner.lexer);
  jq_yylex_destroy(scanner.lexer);
  if (errors > 0) {
    block_free(*answer);
    *answer = gen_noop();
  }
  return errors;
}

int jq_parse_library(struct locfile* locations, block* answer) {
  int errs = jq_parse(locations, answer);
  if (errs) return errs;
  if (!block_has_only_binders(*answer, OP_IS_CALL_PSEUDO)) {
    locfile_locate(locations, UNKNOWN_LOCATION, "error: library should only have function definitions, not a main expression");
    return 1;
  }
  return 0;
}
