function [diffMat, maxDiffMat] = calculateIntegral
	
	clc;
	%clear;

    sRegion = 2;
    sStep = 0.1;
    % sStep = 0.5;
    % pStep = 0.3;
    pStep = 0.1;
    
    S = -sRegion:sStep:sRegion;
    p = 0:pStep:1;

    p(p==0 | p==0.5 | p==1) = [];

	[s_0,s_1,pMesh] = meshgrid(S,S,p);

	newS_0 = s_0(s_0 ~= s_1);
    newS_0 = reshape(newS_0, size(s_0,1), size(s_0,2) - 1, size(s_0,3));

    newS_1 = s_1(s_0 ~= s_1);
    newS_1 = reshape(newS_1, size(s_1,1), size(s_1,2) - 1, size(s_1,3));

    pMesh(:,end,:) = [];

    argDiff = arrayfun(@(a,b,c) myIntegral(a,b,c), newS_0, newS_1, pMesh,'UniformOutput', false);

    % diffMat = struct2table(cell2mat(argDiff(:)));
    diffMat = cell2mat(argDiff(:));

	maxDiffMat = max(diffMat);
end

function result = myIntegral(s_0, s_1, p_0)

	p_1 = 1 - p_0;

	s_pow = (s_0-s_1).^2;

	g_0 = max(p_1,p_0)/min(p_1,p_0);

	g_explicit = @(t) (1 + (1 - 4*p_0*p_1*exp(-s_pow.*t)).^0.5) ./ (1 - (1 - 4*p_0*p_1*exp(-s_pow.*t)).^0.5);
	
	g_integral = @(t) (g_explicit(t) + 1)./(g_explicit(t) - 1);

	% g_tau = @(tau) (arrayfun(@(a,b) integral(g_integral,a,b, 'ArrayValued',true), zeros(size(tau)), tau));
	% g = @(t) g_0.*exp(s_pow.*arrayfun(@(a,b) integral(g_explicit,a,b, 'ArrayValued',true), zeros(size(t)), t));
	g = @(t) g_0.*exp(s_pow.*(arrayfun(@(a,b) integral(g_integral,a,b, 'ArrayValued',true), zeros(size(t)), t)));

	% t = 0:0.1:5;
	t = 0:0.01:0.1;
	% t = 0:0.1:1;
	% t = 0:0.5:1;
    
	g_explicit_result = arrayfun(g_explicit, t);
	g_result = g(t);
    
    result = g_explicit_result - g_result;
end

% t = 0:0.01:0.1;
% 6.66133814775094e-16	5.32907051820075e-15	8.88178419700125e-15	1.06581410364015e-14	7.10542735760100e-15	1.42108547152020e-14	1.24344978758018e-14	2.84217094304040e-14	3.90798504668055e-14	8.52651282912120e-14	4.26325641456060e-14

% t = 0:0.1:5;
% 6.66133814775094e-16	4.26325641456060e-14	1.13686837721616e-12	3.75166564481333e-11	3.84261511499062e-10	2.95112840831280e-08	5.17902662977576e-07	9.02230385690928e-06	3.26978042721748e-05	0.00249296426773071	0.207442179322243	1.47516554594040	21.9320306777954	56.2281427383423	3443.99064636230	482301.968933105	153519.866363525	22297033.7753906	2282766558.77148	47350606990.4297	123741026837.094	157976270581160	935694233352775	689075071837960	133239758422.438	888505985896270	60539672195886.0	532443375528536	256188454359626	1.25758304628511e+15	403193101981.125	2.76805109244465e+15	2.87698774687073e+15	1.14772048980361e+15	3.35371655025656e+15	2.70694356920663e+15	91547303663463.0	712693210563770	49588492345664.0	660456668209417	50933729851876.0	105668403409222	32550007259564.5	481060390689.406	16514453922120.0	263442248397848	311329891069.344	574058986231702	532443375528288	2.48288090287125e+15	16786802507273.0