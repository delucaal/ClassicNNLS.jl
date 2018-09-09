__precompile__()

module ClassicNNLS

export ClassicNNLSSolver

"""
Solves argmin ||X*d-Y||^2 subject to d(i) >= 0.

The original version of this code was developed by
Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
1973 JUN 12, and published in the book
"SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
Revised FEB 1995 to accompany reprinting of the book by SIAM.
"""
function ClassicNNLSSolver(X,Y,MAX_ITER=1000,tol=1e-8)

        N = size(X,2);

        # initialization
        global d = zeros(Float32,N,1);
        global w = transpose(X)*(Y-X*d);
        global R = ones(Bool,N,1);
        global P = zeros(Bool,N,1);
        global mw = maximum(w);

        # Main loop
        for iter in 1:MAX_ITER
            A = findmax(w);
            mw = A[1];
            mIX = A[2];

            # Look for the worst initial solution and include it in the active set
            P[mIX] = true;
            R[mIX] = false;

            # Solve in the active set
            s = zeros(Float32,N,1);
            SEL = P.>0;

            s[SEL] = (X[:,findall(SEL .== true)])\Y;
            # Non-negativity constraint
            neg_subset = s .< 0;
            global d;
            global s;
            global P;
            global R;
            while (sum(neg_subset) > 0)
                neg_subset = (s .< 0);
                alfa = (d[neg_subset]/(d[neg_subset]-s[neg_subset]));
                alfac = findmin(alfa);
                alfac = alfac[1];
                od = d;
                d = d + alfac*(s-d);
                # update the active set
                ToBeUpdated = findall(BitArray(d .== 0.0) .& BitArray(P .== true));

                L2Update = length(ToBeUpdated);
                P[ToBeUpdated] = repeat([false],inner=(1,1),outer=(L2Update,1));
                R[ToBeUpdated] = repeat([true],inner=(1,1),outer=(L2Update,1));

                s = zeros(Float32,N,1);
                SEL = P.>0;
                s[SEL] = (X[:,findall(SEL .== true)])\Y;
                neg_subset = s .< 0;
            end
            d = s;
            w[:] = X'*(Y-X*d);

            if (sum(R) == 0 || mw < tol)
                break
            end
        end

        return d;
end
end #module
