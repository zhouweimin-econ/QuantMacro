% Llnenforce: Enforces LLN in simulating aiyari model

% Here, for a period t the indices have been generated sZ(:,t)
    % Ideally, the number of indices for each status should be equal the
    % ergodic probability. For instance, N=100, pystat=[0.30 0.70], you
    % would need 30 indices of sZ=1 and 70 indices of sZ=2

    
    % Begin with first index, do until the next to last
    % (e.g. with ny=5, do for the first 4 only - 5th is a residual -)
    for j=1:(size(PLR,1)-1)

        % Find the index of all times you have instances of sZ=1
        targets = find(sZ(:,t)==j) ;
        % Find the index of all times you have instances of sZ>1
        % (you don't care about instances smaller than your index,
        % otherwise you would have stealing back and forth)
        misses  = find(sZ(:,t)>j) ;

        % Count the excesses, they are positive is sZ occurs more than 30
        % times, negative viceversa
        excess(j) = round( length(targets) - round(PLR(j)*N) ) ;

        % If the excess are negative, steal from the other indices (the misses) 
        if excess(j)<0
            sZ(misses(1:-excess(j)),t) = j ;
        end

        % If the excesses are positive, give them to the next guy (index j+1)
        if excess(j)>0
            sZ(targets(1:excess(j)),t) = j+1 ;
        end

    end
