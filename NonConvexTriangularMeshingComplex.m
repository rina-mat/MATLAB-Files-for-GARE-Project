function NonConvexTriangularMeshingComplex

format short
tp = 0;
fprintf('Instruction: 1) Enter the nodes as well as the nodal position of the vertices in anticlockwise order \n \t \t 2)The convex subdomains should be entered in anti-clockwise order \n\t\t 3)If number of vertices in a subdomain is less than that of the maximum vertices then fill the rest entries with zeros.\n');
         


n = input('Number of vertices? ');
%------------Input: Domain Co-ordinates------------------------%
fprintf('Enter the co-ordinates: \n');


for i = 1:n
    fprintf('Node %i: ', i);
    node(i, :) = input(''); 
     poly(1,i) = i;
end
lastnode = n+1;



%Concavity/Convexity Check 

[concP] = findConcave(node)
if concP == 0

maxref = input('Enter Number of maximum refinements (Put any integer value (1 to 4 is preferable)): ');

for i = 1:maxref
    ne(i) =  4^(i-1)*n; 
end


oe = input('Enter order of element(1 or 2 or 3):');


 %%%% Functions for meshing elements of different orders of different type %%%

 %%% Linear Element
 
for i = 1: maxref 
if n == 3
   tri = poly;
else
   tri = zeros(1,3);
end
if oe == 1
    
    [nodeCord  triNode] = meshConv1(tri, node, ne(i), n);
    fprintf('For Refinement No.: %i', i);
    print(triNode,nodeCord);
    tri = triNode;
    node = nodeCord;
    nodeCord = node;
    drawDomain(tri, node);
end

%%% Quadratic element of Lagrange Type %%%
if oe == 2 

    [nodeCord  triNode] = meshConv1(tri, node, ne(i), n); 
    tri = triNode;
    node = nodeCord;
   [nodecord  triNode] = mesh2(tri, node);
    fprintf('For Refinement No.: %i', i);
    print(triNode,nodecord);
    drawDomain(triNode,nodecord);
end
%%% Cubic element of Lagrange Type %%%
if oe == 3 
    [nodeCord  triNode] = meshConv1(tri, node, ne(i), n); 
     tri = triNode;
    node = nodeCord;
   [nodecord  triNode] = mesh3(tri, node);
    fprintf('For Refinement No.: %i', i);
    print(triNode,nodecord);
   drawDomain(triNode,nodecord);   
end

end


else

%------------Input: convex subdomain selected by user --------------------%
fprintf('How many parts? \n');
np = input('');
fprintf('Enter the number of maximum vertices: \n');
max = input('');
part = zeros(np,max)
fprintf('Enter each part of the irregular domain: \n');
total = 0;
for i=1:np
    count(i) = 0;
    fprintf('Part[%d]: ',i);

        part(i,:) = input('')
       for j = 1: max
    if part(i,j) ~ 0
            count(i) = count(i)+1;
    end
       end

       count(i)
total = total + count(i);
end


fprintf('Total = %d', total);
display(part);
 
%------------Input: Required number of elements ------------------------%
for i=1:np
    if count(i) < 4
      tp = 1;
      break;
    end
end
ref = input('Enter Number of refinement (Put any value from 1 to 3): ');


sum = 0;
   for j = 1:np
    sum = sum + 4^(ref-1)*count(j);
   end
   
%------------Input: Required order of each elements ---------------------%
oe = input('Enter order of element(1--> Linear or 2--> Quadratic or 3 --> Cubic):');
ne = sum;
cp = np;


%------------Case 1: Order of element is 1 -------------------------%

   
if oe == 1
    tri = zeros(1,3); 
    for i = 1:cp
    if i == 1
        l = 0;
        p = 1;
    end
    [nodeCord  TriNode, ln, pos] = meshConc1(part(i,:), tri, node, ne, cp, count(i), count(i), lastnode, p);
    print(TriNode,nodeCord);
    tri = TriNode;
    node = nodeCord;
    lastnode = ln;
    nodeCord = node;
    p = pos; 
    fprintf('ne = %d, total = %d', ne, total);
    if i == cp && ne ~= total
    fprintf('Accessing meshelm ne = %d, total = %d', ne, total); 
    [nodeCord  TriNode] = meshelm(tri, node, ne, n)
    print(TriNode,nodeCord);
    tri = TriNode;
    node = nodeCord;
    end
    end
    
 drawDomain(tri, node);
end

%------------Case 2: Order of element is 2 -------------------------%
if oe == 2 
tri = zeros(1,3); 
    for i = 1:cp
    if i == 1
        l = 0;
        p = 1;
    end
    [nodeCord  TriNode, ln, pos] = meshConc1(part(i,:), tri, node, ne, cp, count(i), count(i), lastnode, p);
    print(TriNode,nodeCord);
    tri = TriNode;
    node = nodeCord;
    lastnode = ln;
    nodeCord = node;
    p = pos; 
    fprintf('ne = %d, total = %d', ne, total);
   
    %fprintf('i = %d, cp = %d, ne = %d, total = %d', i, cp, ne, total)
    if i == cp && ne == total
    [nodecord  TriNode] = mesh2(tri, node);
    print(TriNode,nodecord);
    drawDomain(TriNode,nodecord);
    end
    if i == cp && ne ~= total
        
    [nodeCord  TriNode] = meshelm(tri, node, ne, n)
   
    
    print(TriNode,nodeCord);
    tri = TriNode;
    node = nodeCord;
    [nodecord  TriNode] = mesh2(tri, node);
    print(TriNode,nodecord);
    drawDomain(TriNode,nodecord);
    end
    end
end




%------------Case 3: Order of element is 3 -------------------------%

if oe == 3
       tri = zeros(1,3); 
    for i = 1:cp
    if i == 1
        l = 0;
        p = 1;
    end
    [nodeCord  TriNode, ln, pos] = meshConc1(part(i,:), tri, node, ne, cp, count(i), count(i), lastnode, p);
    print(TriNode,nodeCord);
    tri = TriNode;
    node = nodeCord;
    lastnode = ln;
    nodeCord = node;
     p = pos; 
    fprintf('ne = %d, total = %d', ne, total)
   
    fprintf('i = %d, cp = %d, ne = %d, total = %d', i, cp, ne, total)
    if i == cp && ne == total
    [nodecord  TriNode] = mesh3(tri, node);
    print(TriNode,nodecord);
    drawDomain(TriNode,nodecord);
    end
    
 if i == cp && ne ~= total
        
    [nodeCord  TriNode] = meshelm(tri, node, ne, n)
   
    %fprintf('before using order 3')
    print(TriNode,nodeCord);
    tri = TriNode;
    node = nodeCord;
    [nodecord  TriNode] = mesh3(tri, node);
    print(TriNode,nodecord);
    drawDomain(TriNode,nodecord);
    end
    end
  


end
end
end

%-----Checking and finding Concave points ------%
    
function [concP] = findConcave(node)

[r c] = size(node)
 q = 0
for i = 1:r  
    if (i == r-1)
        j = r
        k = 1
    elseif (i == r)
            j = 1
            k = 2
        else
            j = i+1
            k = j+1
     end

    
    vx1 = node(j, 1) - node(i,1)
    vy1 = node(j, 2) - node(i,2)
    vx2 = node(k, 1) - node(j,1)
    vy2 = node(k, 2) - node(j,2)
    
    crossProduct(i)= (vx1*vy2 - vx2*vy1)
   
    if crossProduct(i) < 0 
        concave(q+1) = crossProduct(i)
        q = q+1
        p(q) = i+1
    end
end

if (q ~= 0)
fprintf('Polygon Type - Concave \n\n It has concavity at the following points - ')
concP = p
else
    fprintf('Polygon Type - Convex')
    concP = 0;
end

return;
end
 
                       %---Meshing Process: For Case 1------% 
                       
                       
function [node tri] = meshConv1(tri, node, ne, n)
tri
node
ne
if ne == 1
    return;
    
else
    x = [0 0];
    for i = 1:n
        x = x+node(i,:);
    end
    x = x./n
    node(i+1, :) = x
   CentNode = i+1
   lastNode = CentNode 
   
    for j = 1:n
        m = j+1;
        if j==n
            m=1;
        end
         tri(j,:) = [j m CentNode]
    end


if ne == n
    return;
    
else
     for l=1:10
        [noTri g] = size(tri);
        [lastNode g] = size(node);
       
        if noTri == ne
            break;
        end
        elementsFlag = 1;
       % tempNode = zeros(1,4);
        for i=1:noTri
             c = [0 0];
            for j=1:3
                m = j+1;
                if j==3
                    m=1;
                end
                
 x(1) = (node(tri(i,j),1)+node(tri(i,m),1))/2;
 
 x(2) = (node(tri(i,j),2)+node(tri(i,m),2))/2;
 x = [x(1) x(2)];

temp = find(ismember(node,x,'rows'));

if temp<=lastNode
    tempNode(j) = temp;
    c = c + node(temp, :);
    
else
    node(lastNode+1, :) = x;
    lastNode = lastNode + 1;
    tempNode(j) = lastNode;
    c = c + node(lastNode, :);
end

            end
 
    elements(elementsFlag, :) = [tri(i,1) tempNode(1) tempNode(3)];
    elements(elementsFlag+1, :) = [tempNode(1) tri(i,2) tempNode(2)];
    elements(elementsFlag+2, :) = [tempNode(3) tempNode(2) tri(i,3)]; 
    elements(elementsFlag+3, :) = [tempNode(1) tempNode(2) tempNode(3)];
    elementsFlag = elementsFlag + 4;
            
        end    
 tri = elements;

    end

    end
end
end
                       
function [node tri ln pos] = meshConc1(part, tri, node, e, cp, ne, n, lastnode,  p)

pos = p;
if e == cp
    tri(pos,:) = part;
    pos = pos+1;
    ln = lastnode;
    node = node;
    return
end
    x = [0 0]
    fprintf ('ne = %d, n = %d pos = %d', ne, n, pos);
    for i = 1:n
        x = x+node(part(1,i),:);
    end
    x = x./n;
    node(lastnode, :) = x;
   CentNode = lastnode;
   lastnode = CentNode +1;
     for j = 1:n
        m = j+1;
        if j==n
            m=1;
        end
         tri(pos,:) = [part(1,j) part(1,m) CentNode];
         pos = pos +1; 
    end
   

 ln = lastnode;

end

         %---Meshing Process: meshing each subdomain created from the problem domain ------%
         
function [node  tri] =  meshelm(tri, node, ne, n)
    for l=1:10
        [noTri g] = size(tri);
        [lastNode g] = size(node);
       
        if noTri == ne
            break;
        end
        elementsFlag = 1;
       % tempNode = zeros(1,4);
        for i=1:noTri
             c = [0 0];
            for j=1:3
                m = j+1;
                if j==3
                    m=1;
                end
                
 x(1) = (node(tri(i,j),1)+node(tri(i,m),1))/2;
 
 x(2) = (node(tri(i,j),2)+node(tri(i,m),2))/2;
 x = [x(1) x(2)];

temp = find(ismember(node,x,'rows'));

if temp<=lastNode
    tempNode(j) = temp;
    c = c + node(temp, :);
    
else
    node(lastNode+1, :) = x;
    lastNode = lastNode + 1;
    tempNode(j) = lastNode;
    c = c + node(lastNode, :);
end

            end
 
    elements(elementsFlag, :) = [tri(i,1) tempNode(1) tempNode(3)];
   elements(elementsFlag+1, :) = [tempNode(1) tri(i,2) tempNode(2)];
   elements(elementsFlag+2, :) = [tempNode(3) tempNode(2) tri(i,3)]; 
   elements(elementsFlag+3, :) = [tempNode(1) tempNode(2) tempNode(3)];
    elementsFlag = elementsFlag + 4;
            
        end    
 tri = elements;

    end

end

%---Meshing Process: For Case 2------% 

function [node tri] = mesh2(tri,node)
[ne tem] = size(tri);
[l temp2] = size(node);
for e = 1:ne
    temp = (node(tri(e,1), :) + node(tri(e,2),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
    temp = (node(tri(e,2), :) + node(tri(e,3),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end

     temp = (node(tri(e,3), :) + node(tri(e,1),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
     
    
    newTri(e,:) = [tri(e,1) q1 tri(e,2) q2 tri(e,3) q3];
end
tri = newTri;
return;


end



                        %---Meshing Process: For Case 3------% 
                        
                        
function [node tri] = mesh3(tri,node)
[ne tem] = size(tri)
[l temp2] = size(node)
for e = 1:ne
    temp = (2*node(tri(e,1), :) + node(tri(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
     temp = (node(tri(e,1), :) + 2*node(tri(e,2),:))/3;
     q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
    temp = (2 * node(tri(e,2), :) + node(tri(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
 
    temp = (node(tri(e,2), :) + 2*node(tri(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
     temp = (2*node(tri(e,3), :) + node(tri(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
      temp = (node(tri(e,3), :) + 2* node(tri(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
       
    node(l+1, :) = (node(tri(e,1), :) + node(tri(e,2), :) + node(tri(e,3), :))/3; 
    q7 = l+1;
    l = l+1;
    
    newtri(e,:) = [tri(e,1) q1 q2 tri(e,2) q3 q4 tri(e,3) q5 q6 q7];
end
tri = newtri;
return;

end

  %---Output 1: Prints the meshed elements and nodal coordinates of each nodes within these elements ------%


function print(tri, node)
[m n] = size(tri);
[k l] = size(node);
disp('Global Nodes (element wise): ');

for i=1:m
    fprintf('\n[%2i]', i);
    for j = 1:n
        fprintf('%3i ',tri(i,j));
    end
end

fprintf('\nNodal Co-ordinates: \n\n');
fprintf('Node NO.:  abscissa  ordinate\n');

for i=1:k
    fprintf('  %3i : %0.3f    %0.3f\n', i, node(i,1), node(i,2));

end
end


    %---Output 2: Plots the figure of the meshed domain------%

     function drawDomain(t, c)
     [k l] = size(t);
                [r s] = size(c);
                nodesInElm = l;
               % For 3rd order element
                if rem(l,3) ~= 0
                    p = rem(l,3); 
                    l = l - p; 
                end
                
        
        xmin = min(c(:,1)) ;
        ymin = min(c(:,2));
        xmax = max(c(:,1)) ;
        ymax = max(c(:,2));
        
        figure
        xlim([xmin-0.5, xmax+0.5]); ylim([ymin-0.5, ymax+0.5]);  
        hold on ; 
        
        for i=1:k
            cx = sum(c(t(i,:), 1))/ nodesInElm;
            cy = sum(c(t(i,:), 2))/ nodesInElm;
            for j = 1:l
                m = j+1;
                if j == l 
                    m = 1;
                end
                x = [c(t(i,j),1) c(t(i,m),1)];
                y = [c(t(i,j),2) c(t(i,m),2)];
                plot(x,y);
                hold on;  
            end     
          
             %text(cx-0.06,cy,'(');
            % hold on
            text(cx+0.02,cy,int2str(i));
            hold on
            %text(cx-0.02,cy,')');
            %hold on
            
            
        end
          
        for i = 1:k
            for j = 1:l
                plot(c(t(i,j),1), c(t(i,j),2), 'Marker', 'O', 'MarkerFaceColor', 'blue')
                hold on;
            end
            if rem( nodesInElm,3) ~= 0
                  plot(c(t(i, nodesInElm),1), c(t(i, nodesInElm),2), 'Marker', 'O', 'MarkerFaceColor', 'red')
           end
        end
          
          
      title (['Mesh with  ' num2str(k) '  ' num2str(l) '-triangular elements and  ' num2str(r) '  nodes']);
        xlabel('x axis');
        ylabel('y axis');


        
     end