function [c, rec, det, lam, corm] = rqa(gt, sim, r)
    l = min(length(gt), length(sim));

    c = zeros(l, l);

    gt = gt(1:l-1,:);
    sim = sim(1:l-1,:);

    for i = 1:l
        for j = 1:l
            if (sqrt(sum((gt(i) - sim(j)) .^2, 2)) <= r)
                c(i, j) = 1;
            end
        end
    end
    
    rec = crossrec(c);
    det = determinism(c);
    lam = laminarity(c);
    corm = centerrecmass(c);
end

function rec  = crossrec(c)
    C = sum(c(:));
    N = length(c);
    rec = 100*(C/N^2);
end

function det  = determinism(c)
    if any(c(:))
        C = sum(c(:));
        N = length(c);
        D = 0;

        for i = -N:N
            d = sum(diag(c,i));
            if (d > 1)
                D = D + d;
            end
        end

        det = 100*(D/C);
    else
        det = 0;
    end
end

function lam = laminarity(c)
    C = sum(c(:));
    N = length(c); 
    H = 0;
    V = 0;
    
    for i = 1:N
        h = sum(c(i,:));
        v = sum(c(:,i));
        if (h > 1)
            H = H + h;
        end
        if (v > 1)
            V = V + v;
        end
    end
    
    lam = 100*((H+V)/2*C);
end

function corm = centerrecmass(c)
    if any(c(:))
        N = length(c);
        C = sum(c(:));
        Cc = 0;

        for i = 1:N
            for j = 1:N
                Cc = Cc + ((j-i) * c(i,j));
            end
        end

        corm = 100*(Cc/((N-1)*C));
    else
        corm = 0;
    end
end