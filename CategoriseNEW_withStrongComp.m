function output = CategoriseNEW_withStrongComp( Rates_vec, tol, jj)
%function inputs: a vector of Rates
M1fluxE1E2 = Rates_vec(1);   %M1fluxE1E2 = flux of M1 grown in shared environment (no M2)
M1_stM2nw = Rates_vec(2);    %M1fluxE1E2_st_M2NoWorse = flux of M1 grown with M2 growth rate no worse than alone
M1_stM1nw = Rates_vec(3);    %flux of M1 when calculating: M2fluxE1E2_st_M1NoWorse = flux of M2 grown with M1 growth rate no worse than alone
M2fluxE1E2 = Rates_vec(4);   %M2fluxE1E2 = flux of M2 grown in shared environment (no M1)
M2_stM2nw = Rates_vec(5);    %flux of M2 when calculating: M1fluxE1E2_st_M2NoWorse = flux of M1 grown with M2 growth rate no worse than alone
M2_stM1nw = Rates_vec(6);    %M2fluxE1E2_st_M1NoWorse = flux of M2 grown with M1 growth rate no worse than alone
%tol gives the 'close to 0' tolerance 

%output is a scalar value corresponding to the category of the paired
%growth:
%17) No growth (at least one cannot grow)
%%units look at M1 and M2 fluxes when the other does at least no worse:
% (1: +/+), (2: +/0), (3: 0/+), (4: -/-), (5: -/0), (6: 0/-), (7: 0,0)
% (8: +/-), (9: -/+)
%tens: together cases; again, only look at M1|M2nw and M2|M1nw
% (10: T/+), (11: T/0), (12: T/-), (13: +/T), (14:, 0/T), (15: -/T), (16: T/T)
%add 20 to output any time that M1|M1 > M1alone or M2|M2 > M2alone

all_growth = sum( M1fluxE1E2 + M1_stM2nw + M1_stM1nw + M2fluxE1E2 + M2_stM1nw + M2_stM2nw );

if all_growth < tol %no growth
    cat = 17;
elseif M1_stM1nw<tol & M1_stM1nw>-tol & M1_stM2nw<tol & M1_stM2nw>-tol %M1 cannot grow
    cat = 17; 
elseif M2_stM1nw<tol & M2_stM1nw>-tol & M2_stM2nw<tol & M2_stM2nw>-tol %M2 cannot grow
    cat = 17; 
elseif M1fluxE1E2 < tol & M2fluxE1E2 < tol & M1_stM2nw>M1fluxE1E2+tol & M2_stM1nw>M2fluxE1E2+tol %M1 cannot grow alone AND M2 cannot grow alone (but they can together)
        cat = 16;
elseif M1fluxE1E2 < tol %M1 cannot grow alone
        if M2_stM1nw>M2fluxE1E2+tol
            cat = 10;
        elseif M2_stM1nw<=M2fluxE1E2+tol & M2_stM1nw>=M2fluxE1E2-tol
            cat = 11;
        else 
            cat = 12;
        end
elseif M2fluxE1E2 < tol %M2 cannot grow alone
        if M1_stM2nw>M1fluxE1E2+tol
            cat = 13;
        elseif M1_stM2nw<=M1fluxE1E2+tol & M1_stM2nw>=M1fluxE1E2-tol
            cat = 14;
        else 
            cat = 15;
        end
elseif M1_stM2nw>M1fluxE1E2+tol & M2_stM1nw>M2fluxE1E2+tol %both do better when the other is no worse
        cat = 1;
elseif M1_stM2nw>M1fluxE1E2+tol & M2_stM1nw<=M2fluxE1E2+tol & M2_stM1nw>=M2fluxE1E2-tol
        cat = 2;
elseif M1_stM2nw<=M1fluxE1E2+tol & M1_stM2nw>=M1fluxE1E2-tol & M2_stM1nw>M2fluxE1E2+tol
        cat = 3;
elseif M1_stM2nw<M1fluxE1E2-tol & M2_stM1nw<M2fluxE1E2-tol
        cat = 4;
elseif M1_stM2nw<M1fluxE1E2-tol & M2_stM1nw<=M2fluxE1E2+tol & M2_stM1nw>=M2fluxE1E2-tol
        cat = 5;
elseif M1_stM2nw<=M1fluxE1E2+tol & M1_stM2nw>=M1fluxE1E2-tol & M2_stM1nw<M2fluxE1E2-tol
        cat = 6;
elseif M1_stM2nw<=M1fluxE1E2+tol & M1_stM2nw>=M1fluxE1E2-tol & M2_stM1nw<=M2fluxE1E2+tol & M2_stM1nw>=M2fluxE1E2-tol
        cat = 7;
elseif M1_stM2nw>M1fluxE1E2+tol & M2_stM1nw<M2fluxE1E2-tol
        cat = 8;
elseif M1_stM2nw<M1fluxE1E2-tol & M2_stM1nw>M2fluxE1E2+tol
        cat = 9;
else
        str = "something else row " + jj;
        disp(str)
end

output.cat = cat;

%check for strong competition
isStrongComp = 0;
if cat==4
    isStrongComp = M1_stM2nw<tol | M2_stM1nw<tol;
end

output.isStrongComp = isStrongComp;

%     %might be interesting to have a way of documenting whether the
%     %'no worse' growth of a microbe is improved when the other microbe's
%     %growht is being optimised
%     if M1_stM1nw > M1fluxE1E2+tol
%         output = output + 20;
%     end
%     
%     if M2_stM2nw > M2fluxE1E2 + tol
%         output = output + 20;
%     end

end