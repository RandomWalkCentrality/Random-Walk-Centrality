function cholesky()
    fid = fopen('./filename_ic.txt', 'r');
    nn = str2double(fgetl(fid));

    % Create a file to store results
    if exist('decomposition_results.txt', 'file') == 2
        result_file = fopen('decomposition_results.txt', 'a');
    else
        result_file = fopen('decomposition_results.txt', 'w');
        fprintf(result_file, 'Graph Name\tNodes\tEdges\tDecomposition Time (s)\n');
    end

    for i = 1:nn
        str_line = fgetl(fid);
        if contains(str_line, '#')
            continue;
        end

        fprintf('Processing graph: %s\n', str_line);

        G = get_graph(str_line);
        
        % Find the largest connected component
        G_largest = findconnect(G);
        fprintf('Largest connected component found. Nodes: %d, Edges: %d\n', G_largest.n, G_largest.m);

        [norm_L, v] = normalized_laplacian(G_largest);
        norm_L = norm_L(setdiff(1:G_largest.n, v), setdiff(1:G_largest.n, v));
        norm_L = norm_L + speye(size(norm_L)) * 1e-6;

        try
            % Use ichol for incomplete Cholesky decomposition and time it
            opts.type = 'ict';
            opts.droptol = 1e-3;
            tic;
            L = ichol(norm_L, opts);
            decomp_time = toc;

            % Check and replace NaN values
            nan_replaced = false;
            if any(isnan(L(:)))
                L(isnan(L)) = 1;
                nan_replaced = true;
                fprintf('NaN values detected and replaced with 1.\n');
            end

            % Save sparse matrix to file
            save(sprintf('ichol_factor_%s.mat', str_line), 'L', '-v7.3');

            fprintf('Incomplete Cholesky decomposition completed and saved to ichol_factor_%s.mat\n', str_line);
            
            % Write results to file
            fprintf(result_file, '%s\t%d\t%d\t%.6f\n', str_line, G_largest.n, G_largest.m, decomp_time);
        catch e
            fprintf('Failed to compute Incomplete Cholesky decomposition for %s: %s\n', str_line, e.message);
            % If decomposition fails, still record graph info with time marked as N/A
            fprintf(result_file, '%s\t%d\t%d\tN/A\n', str_line, G_largest.n, G_largest.m);
        end

        fprintf('\n');  % Add a blank line between graphs
    end
    fclose(fid);
    fclose(result_file);
end

% findconnect function
function G_largest = findconnect(G)
    n = G.n;
    bcj = zeros(1, G.n);
    idt = 0;
    noc = zeros(1, G.n);
    
    while min(bcj) == 0
        idt = idt + 1;
        label = find(bcj == 0, 1);
        b = zeros(1, n);
        b(1) = label;
        bcj(label) = idt;
        noc(idt) = 1;
        f = 1;
        r = 2;
        while f < r
            for i = 1:length(G.nbr{b(f)})
                if bcj(G.nbr{b(f)}(i)) == 0
                    b(r) = G.nbr{b(f)}(i);
                    bcj(G.nbr{b(f)}(i)) = idt;
                    noc(idt) = noc(idt) + 1;
                    r = r + 1;
                end
            end
            f = f + 1;
        end
    end
    
    [~, cop] = max(noc);
    
    nc = 0;
    Label = containers.Map('KeyType', 'int64', 'ValueType', 'int64');
    Origin = containers.Map('KeyType', 'int64', 'ValueType', 'int64');
    
    function id = getID(x)
        if ~isKey(Label, x)
            nc = nc + 1;
            Label(x) = nc;
        end
        id = Label(x);
    end
    
    mc = 0;
    uc = [];
    vc = [];
    for i = 1:G.m
        if bcj(G.u(i)) == cop
            mc = mc + 1;
            u1 = getID(G.u(i));
            v1 = getID(G.v(i));
            Origin(u1) = G.u(i);
            Origin(v1) = G.v(i);
            uc(end+1) = u1;
            vc(end+1) = v1;
        end
    end
    
    nbr = cell(1, nc);
    for i = 1:mc
        u1 = uc(i);
        v1 = vc(i);
        nbr{u1} = [nbr{u1}, v1];
        nbr{v1} = [nbr{v1}, u1];
    end
    
    degree = cellfun(@length, nbr);
    w = ones(1, mc);
    
    G_largest = struct('n', nc, 'm', mc, 'u', uc, 'v', vc, 'w', w, 'nbr', {nbr}, 'name', G.name, 'degree', degree);
end

function G = get_graph(ffname)
fname = sprintf('./graphs/%s.txt', ffname);
fid = fopen(fname, 'r');

n = 0;
Label = containers.Map('KeyType', 'int64', 'ValueType', 'int64');
Origin = containers.Map('KeyType', 'int64', 'ValueType', 'int64');

u = [];
v = [];

while ~feof(fid)
str = fgetl(fid);
if isempty(str) || str(1) == '#'
continue;
end

if contains(str, ',')
    data = sscanf(str, '%d,%d');
else
    data = sscanf(str, '%d %d');
end

if length(data) ~= 2
    continue;
end

x = data(1);
y = data(2);

if ~isKey(Label, x)
    n = n + 1;
    Label(x) = n;
    Origin(n) = x;
end
if ~isKey(Label, y)
    n = n + 1;
    Label(y) = n;
    Origin(n) = y;
end

u(end+1) = Label(x);
v(end+1) = Label(y);
end

fclose(fid);

m = length(u);
nbr = cell(1, n);
for i = 1:m
nbr{u(i)} = [nbr{u(i)}, v(i)];
nbr{v(i)} = [nbr{v(i)}, u(i)];
end

degree = cellfun(@length, nbr);
w = ones(1, m);

G = struct('n', n, 'm', m, 'u', u, 'v', v, 'w', w, 'nbr', {nbr}, 'name', ffname, 'degree', degree);
end

function [norm_L, max_degree_node] = normalized_laplacian(G)
n = G.n;
rows = [G.u, G.v, 1:n];
cols = [G.v, G.u, 1:n];
data = [-ones(1, G.m), -ones(1, G.m), G.degree];

fprintf('Length of rows: %d\n', length(rows));
fprintf('Length of cols: %d\n', length(cols));
fprintf('Length of data: %d\n', length(data));
fprintf('n: %d, m: %d\n', n, G.m);

if ~(length(rows) == length(cols) && length(cols) == length(data))
error('Arrays must have the same length. Rows: %d, Cols: %d, Data: %d', length(rows), length(cols), length(data));
end

L = sparse(rows, cols, data, n, n);

zero_degree = (G.degree == 0);
deg_sqrt_inv = G.degree.^(-0.5);
deg_sqrt_inv(zero_degree) = 0;

deg_sqrt_inv = deg_sqrt_inv(:);

D_sqrt_inv = spdiags(deg_sqrt_inv, 0, n, n);

norm_L = D_sqrt_inv * L * D_sqrt_inv;

[~, max_degree_node] = max(G.degree);
end