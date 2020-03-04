function [all_Mat, diff_Mat] = calculateIntegral
	
	clc;
	clear;

	sRegion = 2;
	
	sStep = 0.1;
	% sStep = 0.5;
	pStep = 0.1;
	% pStep = 0.3;
	
	S = -sRegion:sStep:sRegion;
	p = 0:pStep:1;

	p(p==0 | p==0.5 | p==1) = [];

	[s_0,s_1,pMesh] = meshgrid(S,S,p);

	newS_0 = s_0(s_0 ~= s_1);
	newS_0 = reshape(newS_0, size(s_0,1), size(s_0,2) - 1, size(s_0,3));

	newS_1 = s_1(s_0 ~= s_1);
	newS_1 = reshape(newS_1, size(s_1,1), size(s_1,2) - 1, size(s_1,3));

	pMesh(:,end,:) = [];

	[all_ret, all_diff] = arrayfun(@(a,b,c) myIntegral(a,b,c), newS_0, newS_1, pMesh,'UniformOutput', false);

	all_Mat = cell2mat(all_ret(:));
	diff_Mat = cell2mat(all_diff(:));
end

function [all, all_diff] = myIntegral(s_0, s_1, p_0)

	p_1 = 1 - p_0;

	s_pow = (s_0-s_1).^2;

	g_0 = max(p_1,p_0)/min(p_1,p_0);

	g_explicit = @(t) (1 + (1 - 4*p_0*p_1*exp(-s_pow.*t)).^0.5) ./ (1 - (1 - 4*p_0*p_1*exp(-s_pow.*t)).^0.5);
	
	g_integral = @(t) (g_explicit(t) + 1)./(g_explicit(t) - 1);

	g = @(t) g_0.*exp(s_pow.*(arrayfun(@(a,b) integral(g_integral,a,b, 'ArrayValued',true), zeros(size(t)), t)));

	t = 0:0.1:5;
	% t = 0:0.1:1;
	% t = 0:0.01:0.1;
	
	g_explicit_result = arrayfun(g_explicit, t);
	g_result = g(t);
	
	argDiff = g_explicit_result - g_result;
	
	all.diff = argDiff;
	all_diff = argDiff;
	all.s_0 = s_0;
	all.s_1 = s_1;
	all.p_0 = p_0;
end