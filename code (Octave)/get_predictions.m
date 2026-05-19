function [Res] = get_predictions( params, R, SIR, SIC)
    
    gamma = params(1);
    bethas = params(2:1+length(SIR));
    alphas = params(2+length(SIR):length(params));

    Res = our_method(R, gamma, SIR, SIC , alphas, bethas);

end