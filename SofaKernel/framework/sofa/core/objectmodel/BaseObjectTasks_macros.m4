define(`M4_TMP1',`ifelse($1,1,define(`M4_TMP',`$2')M4_TMP($1),`M4_TMP1(decr($1),`$2',`$3')'`$3'define(`M4_TMP',`$2')M4_TMP($1))')
define(`M4_PARAM',`ifelse(KAAPI_NUMBER_PARAMS,0,,`$2M4_TMP1(KAAPI_NUMBER_PARAMS,``$1'',``$3'')$4')')
define(`M4_PARAMV',`ifelse($4,0,,`$2M4_TMP1($4,``$1'',``$3'')$5')')
define(`KAAPI_CLOSUREFMT',ClosureFormat$1)
define(`KAAPI_CLOSURE',KaapiClosure$1)
define(`KAAPI_VERIF',VerifParamType$1)
define(`KAAPI_INITFORMATCLOSURE',InitFormatClosure$1)
define(M4_PARAM_0,`ifelse(KAAPI_NUMBER_PARAMS,0,$@,)')
