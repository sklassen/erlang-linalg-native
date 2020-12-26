-module(linalg_arithmetric_tests). 
-import(linalg_arithmetric,[exp/1,log/1,norm/1]).
-import(linalg_arithmetric,[mul/2,add/2]).
-include_lib("eunit/include/eunit.hrl").

exp_s_test() ->
	?assert(exp(1)==math:exp(1)).
exp_v_test() ->
	?assert(exp([1,2,3])==[math:exp(X)||X<-[1,2,3]]).
exp_m_test() ->
	?assert(exp([[1,2,3],[1,2,3]])==[[math:exp(X)||X<-[1,2,3]]||_<-[1,2]]).

log_s_test() ->
	?assert(log(1)==math:log(1)).
log_v_test() ->
	?assert(log([1,2,3])==[math:log(X)||X<-[1,2,3]]).
log_m_test() ->
	?assert(log([[1,2,3],[1,2,3]])==[[math:log(X)||X<-[1,2,3]]||_<-[1,2]]).

norm_s_test() ->
	?assertEqual(1,norm(1)).
norm_v_test() ->
	?assertEqual(math:sqrt(14),norm([1,2,3])).
norm_m_test() ->
	?assertEqual(math:sqrt(28),norm([[1,2,3],[1,2,3]])).

mul_sxs_test() ->
	?assert(mul(10,10)==100).
mul_sxv_test() ->
	?assert(mul(10,[1,2,3])==[10,20,30]).
mul_vxs_test() ->
	?assert(mul([4,5,6],10)==[40,50,60]).
mul_vxv_test() ->
	?assert(mul([1,2,3],[4,5,6])==[4,10,18]).
mul_sxm_test() ->
	?assert(mul(10,[[1,2,3],[4,5,6]])==[[10,20,30],[40,50,60]]).
mul_mxs_test() ->
	?assert(mul([[1,2,3],[4,5,6]],10)==[[10,20,30],[40,50,60]]).
mul_1x1_test() ->
	?assert(mul([[3]],[[3]])==[[9]]).
mul_mxm_test() ->
	?assert(mul([[1,2,3],[4,5,6]],[[1,2,3],[4,5,6]])==[[1,4,9],[16,25,36]]).

add_test() ->
    [
    ?assert(add(2, 3) == 5),
    ?assert(add([1, 2], [3, 4]) == [4, 6]),
    ?assert(add([[1, 2], [3, 4]], [[2, 3], [4, 5]]) == [[3, 5], [7, 9]]),
    ?assert(add(10, [[2, 3], [4, 5]]) == [[12, 13], [14, 15]])
    ].
