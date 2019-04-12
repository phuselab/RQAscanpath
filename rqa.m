function [c, rec, det, lam, corm] = rqa(s1, s2, r, ll)
% Parameters:
%
% s1                    (n,2) matrix of x,y fixations
% s2                    (n,2) matrix of x,y fixations
% r                     e.g. 34, radius threshold
% ll                    e.g. 2, minimum line length
%
% Results:
%
% c                     Recurrence matrix
% rec                   Cross-recurrence
% det                   Determinism
% lam                   Laminarity
% corm                  Center of recurrence mass
%
% Author                Vittorio Cuculo
% Date                  April 2019

% Determine the shorter length
l = min(length(s1), length(s2));

% Initialize recurrence matrix
c = zeros(l, l);

% Truncate the sequences to the shortest
s1 = s1(1:l-1,:);
s2 = s2(1:l-1,:);

% Compute recurrence matrix
for i = 1:l
    for j = 1:l
        if (pdist2(s1(i),s2(j)) <= r)
            c(i, j) = 1;
        end
    end
end

if any(c(:))
    rec = crossrec(c);
    det = determinism(c, ll);
    lam = laminarity(c, ll);
    corm = centerrecmass(c);
else
    rec = NaN;
    det = NaN;
    lam = NaN;
    corm = NaN;
end
end

function rec = crossrec(c)
% Sum of recurrences
C = sum(c(:));
% Number of fixations
N = length(c);
rec = 100*(C/N^2);
end

function det  = determinism(c, ll)
C = sum(c(:));
N = length(c);
D = 0;
for i = -N:N
    % Cardinality of diagonal line
    d = sum(diag(c,i));
    % Consider if long at least "line length"
    if (d >= ll)
        D = D + d;
    end
end
det = 100*(D/C);
end

function lam = laminarity(c, ll)
C = sum(c(:));
N = length(c);
H = 0;
V = 0;
for i = 1:N
    % Cardinality of horizontal line
    h = sum(c(i,:));
    % Cardinality of vertical line
    v = sum(c(:,i));
    % Consider if long at least "line length"
    if (h >= ll)
        H = H + h;
    end
    if (v >= ll)
        V = V + v;
    end
end
lam = 100*((H+V)/2*C);
end

function corm = centerrecmass(c)
N = length(c);
C = sum(c(:));
Cc = 0;
for i = 1:N
    for j = 1:N
        Cc = Cc + ((j-i) * c(i,j));
    end
end
corm = 100*(Cc/((N-1)*C));
end