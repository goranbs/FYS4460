function [r_new,direction] = walk(z,lx,ly,r)
%
% USEAGE:
%
% [r_new,direction] = walk(z,lx,ly,lz,r_old)
%
% z = rand(lx,ly) < p, for some probability p in range [0,1].
% r = [x,y], some position on z different from 0.
%
% The random walker returns the direction [north,south,east,vest] it walked
% direction = 1,2,3,4
% 1 == south    r(1) += 1
% 2 == north    r(1) -= 1
% 3 == vest     r(2) -= 1
% 4 == east     r(2) += 1
%
% if we meet the wall (r(0,:), r(:,0), r(lx+1,:), r(:,ly), the randomwalker
% choses another direction to walk.
% if we meet obsticle (z(rx,ry) = 0), the randomwalker choses another
% direction to walk.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose a new way
    function [r_new] = Way(r,way)
        % walking n,s,e,v depending on the; way
        % remember that matlab counts from 1,1 in the top north-vest corner
        if (way == 1)
            r(1) = r(1) + 1;  % south
        elseif (way == 2)
            r(1) = r(1) - 1;  % north
        elseif (way == 3)
            r(2) = r(2) - 1;  % vest
        else
            r(2) = r(2) + 1;  % east
        end
        r_new = r;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hitting the wall
    function [val,r] = have_we_met_the_wall(r)
        % Returns val = 0 if we have not met the wall
        % -//-    val = 1 if we have met the wall in north dir.
        % -//-    val = 2 -//-                       south dir.
        % -//-    val = 3 -//-                       vest dir.
        % -//-    val = 4 -//-                       east dir.
        
        
        val = 0;
        % forward / backward
        if ( r(1) == 0 )
            % must walk south!
            r(1) = r(1) + 1;
            val = 1;
        elseif (r(1) == lx+1)
            % must walk north!
            r(1) = r(1) - 1;
            val = 2;
        end
        
        % left / right
        if ( r(2) == 0 )
            % must walk East!
            r(2) = r(2) + 1;
            val = 3;
        elseif (r(2) == ly+1)
            % must walk Vest!
            r(2) = r(2) - 1;
            val = 4;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if we hit an obsticle, we have to choose a different way!
% and wa have to check if we've hit the wall again

    function [r,val] = obsticle(r_new,r_old,z,way)
        % if we meet an obsticle, we have to choose a different way!
        val = 0;
        if (z(r_new(1),r_new(2)) == 0)
            val = way;
            r = r_old; % set position back to original
        end
        r = r_new;     % keep new position
            
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finding the next step r

dir = [1,2,3,4];
way = (1 + round(3*rand(1)));   % 1,2,3,4

r_old = r;
r_new = Way(r,way);
% now we have found a new position, r_new. Is this position allowed?

[val,r_new] = have_we_met_the_wall(r_new); % if val = 0, have not met the wall
dir(dir==val) = [];
[r_new,val2] = obsticle(r_new,r_old,z,way); % if val2 = 0, have not met obsticle
dir(dir==val2) = [];

allowed = 'false';
if (val == 0 && val2 == 0)
    allowed = 'true';
else
    % 1)
    % try new way:
    index = (1 + round((length(dir) - 1)*rand(1))); % random index in dir
    way = dir(index); 
    r_new = Way(r,way);
    
    [val,r_new] = have_we_met_the_wall(r_new); % if val = 0, have not met the wall
    dir(dir==val) = [];
    [r_new,val2] = obsticle(r_new,r_old,z,way); % if val2 = 0, have not met obsticle
    dir(dir==val2) = [];
    
    if (val == 0 && val2 == 0)
        allowed = 'true';
    else
        % 2)
        % try new way:
        index = (1 + round((length(dir) - 1)*rand(1))); % random index in dir
        way = dir(index);
        r_new = Way(r,way);
        
        [val,r_new] = have_we_met_the_wall(r_new); % if val = 0, have not met the wall
        dir(dir==val) = [];
        [r_new,val2] = obsticle(r_new,r_old,z,way); % if val2 = 0, have not met obsticle
        dir(dir==val2) = [];
        
        if (val == 0 && val2 == 0)
            allowed = 'true';
        else
            % 3)
            % try new way:
            index = (1 + round((length(dir) - 1)*rand(1))); % random index in dir
            way = dir(index);
            r_new = Way(r,way);
            
            [val,r_new] = have_we_met_the_wall(r_new); % if val = 0, have not met the wall
            dir(dir==val) = [];
            [r_new,val2] = obsticle(r_new,r_old,z,way); % if val2 = 0, have not met obsticle
            dir(dir==val2) = [];
            
            if (val == 0 && val2 == 0)
                allowed = 'true';
            else
                % all four alternatives must have been tried!
                %Warning = ' somethings not right'
                %Retures = 'Old position returned'
                r_new = r_old;
                allowed = 'false';
            end
        end
    end
    
    
end

% what way did we move?
Rdiff = r_new - r_old;
if (Rdiff(1) == 0)
    % did not move in north / south dir.
    % direction = vest or east
    if (Rdiff(2) == -1)
        direction = 3; % vest
        Dir = 'vest';
    elseif (Rdiff(2) == 1)
        direction = 4; % east
        Dir = 'east';
    else
        % Rdiff(2) = 0
        % did not move
        direction = 0;
        Dir = 'stayed';
    end
elseif (Rdiff(1) == -1)
    direction = 2; % north
    Dir = 'north';
else
    direction = 1; % south
    Dir = 'south';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Need to go through this ->
% while (allowed == false)
%     % 1) have we met the wall?
%     [val,r_new] = have_we_met_the_wall(r);
%     dir(dir==val) = []; % remove possible direction to move in.
%     
%     % 2) have we met an obsticle?
%     [r_new,val2] = obsticle(r_new,r_old,z,way);
%     dir(dir==val2) = []; % remove possible direction to move in.
%     
%     if (length(dir) == 0)
%         Warning = 'Error! No direction to walk!'
%     end
%     
%     % allowed to move on should be false as long as r_old == r_new!
%     
%     if (isequal(r_old,r_new) == false)
%         allowed = true; % we have a new allowed position to walk to
%     else
%         index = (1 + round((length(dir) - 1)*rand(1))); % random index in dir
%         way = dir(index);
%         r_new = Way(r_old,way);
%         %allowed = false;
%     end
% end
% 


end