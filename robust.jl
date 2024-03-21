#using CPLEX
#using JuMP
const n=5000
# le programme nominal est Ax<=1000_m, oÃš A est une 
# matrice n-par-m, est 1000_m est un vecteur colonne 
# avec la valeur 1000 rÃŠpÃŠtÃŠe m fois. La fonction
# objectif est la somme sum(x[i] for i in 1:n).
function loadData(n::Int64, m::Int64)
    A = zeros(m,n);
    A[1,1] = (m+n)%159;
    for i in 2:n
        A[1,i]=(A[1,1]*A[1,i-1]) % 159 ;
    end
    for j in 2:m
        A[j,1] = A[j-1,1]+3;
        for i in 2:n
            A[j,i]=(A[j-1,i]*A[j,i-1])%159;
        end
    end
    return A;
end
#println(loadData(10, 4))