function [i_sort] = match_points(C)
% use linear programming to find optimal assignment of elements in C2 to
% elements in C1

% C = C;

% length of set
n = size(C,1);

% starting point
x0 = eye(n); x0 = x0(:);

%% Form cost vector
% ======================================================================= %
% C = zeros(n,1);
% for i = 1:n
%     C(((i-1)*n+1):i*n) = sum((P2-ones(n,1)*P1(i,:)).^2,2);
% end

C = C(:);

%% single assigment constraints
% ======================================================================= %
% A1eq = zeros(n,n^2);
% A2eq = zeros(n,n^2);
% A1eq = spalloc(n,n^2,n^2);
% A2eq = spalloc(n,n^2,n^2);

% n = 10;

blocks = ones(n,1)*(1:n);
blocks2 = blocks';

Aeq_rows = [blocks(:);blocks2(:)+n];
Aeq_cols = [(1:n^2)';(1:n^2)'];
Aeq_vals = ones(2*n^2,1);
Aeq = sparse(Aeq_rows,Aeq_cols,Aeq_vals);

beq = ones(2*n,1);

% options = optimoptions('linprog','Algorithm','active-set');
% x = linprog(C,[],[],Aeq,beq,0,1,x0,options);
options = optimoptions('linprog','Display','off');
x = linprog(C,[],[],Aeq,beq,zeros(size(x0)),ones(size(x0)),[],options);


%% turn x into sorting index
% ======================================================================= %

[~,i_sort] = max((reshape(x,[n,n])));
