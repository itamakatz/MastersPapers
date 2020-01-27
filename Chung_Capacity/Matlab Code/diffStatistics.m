function [diffMat,M,I] = diffStatistics
    
    lRegion = 11;
    lStep = 0.001;

    sRegion = 2;
    sStep = 0.1;
    
    pStep = 0.01;
    
    L = -lRegion:lStep:lRegion;
    S = -sRegion:sStep:sRegion;
    % S = 0.1:0.1:2;
    p = 0:pStep:1;

    p(p==0 | p==0.5 | p==1) = [];

    % p = 0.28;
    % p = 0.1;

    [s_0,s_1,pMesh] = meshgrid(S,S,p);
    
    newS_0 = s_0(s_0 ~= s_1);
    newS_0 = reshape(newS_0, size(s_0,1), size(s_0,2) - 1, size(s_0,3));

    newS_1 = s_1(s_0 ~= s_1);
    newS_1 = reshape(newS_1, size(s_1,1), size(s_1,2) - 1, size(s_1,3));

    pMesh(:,end,:) = [];

    argDiff = arrayfun(@(a,b,c) argMaxDiff(a,b,c,L), newS_0, newS_1, pMesh,'UniformOutput', false);
    % argDiff = arrayfun(@(a,b,c) argMaxDiff(a,b,c,L), s_0, s_1, pMesh,'UniformOutput', false);
    
    diffMat = struct2table(cell2mat(argDiff(:)));
    
    % diffMat(diffMat.s_0 == diffMat.s_1,:) = [];
    
    [M,I] = max(diffMat.diff) 

    % sizeS = size(S,2);
    % s_0_Reshape = reshape(diffMat.diff, sizeS, sizeS - 1);
    % s_1_Reshape = reshape(diffMat.s_0, sizeS, sizeS - 1);
    % diffReshape = reshape(diffMat.s_1, sizeS, sizeS - 1);
    % surf(s_0_Reshape, s_1_Reshape, diffReshape);

end

function all = argMaxDiff(s_0, s_1, p_0, L)
    
    p_1 = 1 - p_0;
    L_star = (s_0 .* p_0 - s_1 .* p_1) ./ (p_1 - p_0);

    lambda_0 = (s_0 + L).^2;
    lambda_1 = (s_1 + L).^2;

    lambda_avg = lambda_0 .* p_0 + lambda_1 .* p_1;
    lambda_wave = sqrt(lambda_0) .* p_0 + sqrt(lambda_1) .* p_1;

    I = - lambda_avg .* log(lambda_avg) ...
        + p_0 .* lambda_0 .* log(lambda_0) ...
        + p_1 .* lambda_1 .* log(lambda_1);
    
    % I(abs(I) < 1e-15) = 0;

    [maxValue, maxIndex] = max(I);
    argDiff = L(maxIndex) - L_star;

    all.diff = argDiff;
    all.L_max = L(maxIndex);
    all.L_star = L_star;
    all.s_0 = s_0;
    all.s_1 = s_1;
    all.p_0 = p_0;
end


function print(str, val)
    disp([str ':  ' num2str(val)]);
end