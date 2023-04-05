[1, 2, 3, 4][[true, false, false, true]]

a = collect(1:10)

A = [1 2 3
     4 5 6
     7 8 9
     7 8 9]

B = [1 2 ; 3 4]

A[1,:]
A[:,1]

size(A, 1)
size(A, 2)

transpose(A)
A'

filter(a->a>5, a)

x, ___ = size([2 2; 1 1])

b = ("s", 3)

