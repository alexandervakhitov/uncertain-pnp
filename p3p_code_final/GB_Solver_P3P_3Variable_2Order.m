% Generated using GBSolver generator Copyright Martin Bujnak,
% Zuzana Kukelova, Tomas Pajdla CTU Prague 2008.
% 
% Please refer to the following paper, when using this code :
%     Kukelova Z., Bujnak M., Pajdla T., Automatic Generator of Minimal Problem Solvers,
%     ECCV 2008, Marseille, France, October 12-18, 2008

function [x y z] = GB_Solver_P3P_3Variable_2Order(a1, a2, b1, b2, c1, c2, d1, d2, d3)

	% precalculate polynomial equations coefficients
	c(1) = 1+a2^2+a1^2;
	c(2) = -2*a2*b2-2-2*a1*b1;
	c(3) = b2^2+1+b1^2;
	c(4) = -d1;
	c(5) = 1+a2^2+a1^2;
	c(6) = -2*c2*a2-2-2*c1*a1;
	c(7) = c2^2+1+c1^2;
	c(8) = -d2;
	c(9) = 1+b1^2+b2^2;
	c(10) = -2*c2*b2-2-2*c1*b1;
	c(11) = c2^2+1+c1^2;
	c(12) = -d3;

	M = zeros(26, 34);
	ci = [17, 42, 119, 144, 221, 370, 395, 472, 625];
	M(ci) = c(1);

	ci = [43, 68, 145, 170, 247, 396, 421, 498, 651];
	M(ci) = c(2);

	ci = [69, 94, 171, 196, 273, 422, 447, 524, 677];
	M(ci) = c(3);

	ci = [667, 692, 717, 742, 767, 786, 811, 836, 859];
	M(ci) = c(4);

	ci = [22, 47, 124, 149, 226, 373, 398, 475, 626];
	M(ci) = c(5);

	ci = [152, 177, 228, 253, 304, 477, 502, 553, 704];
	M(ci) = c(6);

	ci = [256, 281, 306, 331, 356, 555, 580, 605, 756];
	M(ci) = c(7);

	ci = [672, 697, 722, 747, 772, 789, 814, 839, 860];
	M(ci) = c(8);

	ci = [104, 181, 206, 283, 428, 453, 530, 679];
	M(ci) = c(9);

	ci = [208, 259, 284, 335, 506, 531, 582, 731];
	M(ci) = c(10);

	ci = [286, 311, 336, 361, 558, 583, 608, 757];
	M(ci) = c(11);

	ci = [702, 727, 752, 777, 792, 817, 842, 861];
	M(ci) = c(12);


	Mr = mex_rref(M);  % replace me with a MEX

	A = zeros(8);
	amcols = [34 33 32 31 30 29 28 24];
	A(1, 4) = 1;
	A(2, 7) = 1;
	A(3, :) = -Mr(25, amcols);
	A(4, :) = -Mr(24, amcols);
	A(5, :) = -Mr(22, amcols);
	A(6, :) = -Mr(20, amcols);
	A(7, :) = -Mr(19, amcols);
	A(8, :) = -Mr(12, amcols);

	[V D] = eig(A);
	sol =  V([4, 3, 2],:)./(ones(3, 1)*V(1,:));

	if (find(isnan(sol(:))) > 0)
		
		x = [];
		y = [];
		z = [];
	else
		
		I = find(not(imag( sol(1,:) )));
		x = sol(1,I);
		y = sol(2,I);
		z = sol(3,I);
	end
end
