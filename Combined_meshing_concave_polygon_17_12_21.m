function QuadTriangularMeshing_03_01

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
    ne(i) =  (4^i)*n 
end
oe = input('Enter order of element(1 or 2 or 3):');

%type = input('Type of Element(LAGRANGE = 1, SERENDIPITY = 2)    ');
 %%%% Functions for meshing elements of different orders of different type %%%

 %%% Linear Element
 for i = 1: maxref 
  if n == 4
    quad = poly;
else 
    quad = zeros(1,4);
  end
if oe == 1
    [nodeCord  QuadTri] = meshConv1(quad, node, ne(i), n);
    print(QuadTri,nodeCord);
    quadtri = QuadTri;
    node = nodeCord;
  %  nodeCord = node;
   drawDomain(quadtri, node);
    [rq cq] = size(QuadTri);
    [rn cn] = size(node);


end

if oe == 2 && type == 1
    
    [nodeCord  QuadNode] = meshConv1(quad, node, ne(i), n); 
    quad = QuadNode;
    node = nodeCord;
   [nodecord  QuadNode] = mesh2L(quad, node);
    print(QuadNode,nodecord);
    drawDomain(QuadNode,nodecord);
    [rq cq] = size(QuadNode);
    [rn cn] = size(nodecord);
    
end
if oe == 2 && type == 2
    
    [nodeCord  QuadNode] = meshConv1(quad, node, ne(i), n); 
    quad = QuadNode;
    node = nodeCord;
   [nodecord  QuadNode] = mesh2S(quad, node);
    print(QuadNode,nodecord);
    drawDomain(QuadNode,nodecord);
    [rq cq] = size(QuadNode);
    [rn cn] = size(nodecord); 


end

if oe == 3 && type == 1
    [nodeCord  QuadNode] = meshConv1(quad, node, ne(i), n); 
    quad = QuadNode;
    node = nodeCord;
   [nodecord  QuadNode] = mesh3L(quad, node);
    print(QuadNode,nodecord);
   drawDomain(QuadNode,nodecord);
    [rq cq] = size(QuadNode);
    [rn cn] = size(nodecord);
    
end
if oe == 3 && type == 2
    
      [nodeCord  QuadNode] = meshConv1(quad, node, ne(i), n); 
    quad = QuadNode;
    node = nodeCord;
   [nodecord  QuadNode] = mesh3S(quad, node);
    print(QuadNode,nodecord);
   drawDomain(QuadNode,nodecord);
    [rq cq] = size(QuadNode);
    [rn cn] = size(nodecord);  


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
meshType = input('Mesh type of each part: ');


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


ne = 0;
   for j = 1:np
    ne = ne + 4^(ref-1)*count(j);
   end

%------------Input: Required order of each elements ---------------------%
oe = input('Enter order of element(1 or 2 or 3):');
elm = ne;
cp = np;

type = input('Type of Element(LAGRANGE = 1, SERENDIPITY = 2)    ');
%------------Case 1: Order of element is 1 -------------------------%
if oe == 1
    quad = zeros(1,4);
    for i = 1:cp
    if i == 1
        l = 0;
    end
    [nodeCord  QuadNode, ln, last] = meshConc1(part(i,:), quad, node, elm, cp, count(i), count(i), lastnode, l);
    print(QuadNode,nodeCord);
    quad = QuadNode;
    node = nodeCord;
    lastnode = ln;
    nodeCord = node;
    l = last;
    [rq cq] = size(QuadNode);
    [rn cn] = size(node);
    fprintf('ne = %d, total = %d', ne, total);
    if i == cp && ne ~ total
        
     [nodeCord  QuadNode] = meshelm(quad, node, ne, n)
    print(QuadNode,nodeCord);
    quad = QuadNode;
    node = nodeCord;
    lastnode = ln;
    nodeCord = node;
    
    [rq cq] = size(QuadNode);
    [rn cn] = size(node);
    end
    end
    
 drawDomain(quad, node);
end

%------------Case 2: Order of element is 2 -------------------------%
if oe == 2 && type == 1
     quad = zeros(1,4);
    for i = 1:cp
    if i == 1
        l = 0
    end
    [nodeCord  QuadNode, ln, last] = meshConc1(part(i,:), quad, node, elm, cp, count(i), count(i), lastnode, l);
    print(QuadNode,nodeCord);
    quad = QuadNode;
    node = nodeCord;
    lastnode = ln;
    nodeCord = node;
    l =last;
    [rq cq] = size(QuadNode);
    [rn cn] = size(node);
   
    fprintf('i = %d, cp = %d, ne = %d, total = %d', i, cp, ne, total)
    if i == cp && ne == total
    [nodecord  QuadNode] = mesh2L(quad, node);
    print(QuadNode,nodecord);
    drawDomain(QuadNode,nodecord);
    [rq cq] = size(QuadNode);
    [rn cn] = size(nodecord);
    end
    if i == cp && ne ~ total
        
    [nodeCord  QuadNode] = meshelm(quad, node, ne, n)
   
    
    print(QuadNode,nodeCord);
    quad = QuadNode;
    node = nodeCord;
    [nodecord  QuadNode] = mesh2L(quad, node);
    print(QuadNode,nodecord);
    drawDomain(QuadNode,nodecord);
    [rq cq] = size(QuadNode);
    [rn cn] = size(nodecord);
    end
    end
end


if oe == 2 && type == 2
     quad = zeros(1,4);
    for i = 1:cp
    if i == 1
        l = 0
    end
    [nodeCord  QuadNode, ln, last] = meshConc1(part(i,:), quad, node, elm, cp, count(i), count(i), lastnode, l);
    print(QuadNode,nodeCord);
    quad = QuadNode;
    node = nodeCord;
    lastnode = ln;
    nodeCord = node;
    l =last;
    [rq cq] = size(QuadNode);
    [rn cn] = size(node);
   
    fprintf('i = %d, cp = %d, ne = %d, total = %d', i, cp, ne, total)
    if i == cp && ne == total
    [nodecord  QuadNode] = mesh2S(quad, node);
    print(QuadNode,nodecord);
    drawDomain(QuadNode,nodecord);
    [rq cq] = size(QuadNode);
    [rn cn] = size(nodecord);
    end
    if i == cp && ne ~ total
        
    [nodeCord  QuadNode] = meshelm(quad, node, ne, n)
   
    %fprintf('before using order 2')
    print(QuadNode,nodeCord);
    drawDomain(QuadNode,nodeCord);
    quad = QuadNode;
    node = nodeCord;
    
   %lastnode = ln;
   % nodeCord = node;
    [nodecord  QuadNode] = mesh2S(quad, node);
    print(QuadNode,nodecord);
    drawDomain(QuadNode,nodecord);
    [rq cq] = size(QuadNode);
    [rn cn] = size(nodecord);
    end
    end
  


end

%------------Case 3: Order of element is 2 -------------------------%

if oe == 3 && type == 1
    quad = zeros(1,4);
    for i = 1:cp
    if i == 1
        l = 0;
    end
    [nodeCord  QuadNode, ln, last] = meshConc1(part(i,:), quad, node, elm, cp, count(i), count(i), lastnode, l);
    print(QuadNode,nodeCord);
    quad = QuadNode;
    node = nodeCord;
    lastnode = ln;
    nodeCord = node;
    l = last
    [rq cq] = size(QuadNode);
    [rn cn] = size(node);
    fprintf('ne = %d, total = %d', ne, total);
    if i == cp && ne == total
    [nodecord  QuadNode] = mesh3L(quad, node);
    print(QuadNode,nodecord);
    drawDomain(QuadNode,nodecord);
    [rq cq] = size(QuadNode);
    [rn cn] = size(node);
    end
    
    if i == cp && ne ~ total
        
    [nodeCord  QuadNode] = meshelm(quad, node, ne, n);
    print(QuadNode,nodeCord);
    quad = QuadNode;
    node = nodeCord;
    lastnode = ln;
    [nodecord  QuadNode] = mesh3L(quad, node);
    print(QuadNode,nodecord);
    drawDomain(QuadNode, nodecord);
    [rq cq] = size(QuadNode);
    [rn cn] = size(nodecord);
    end
    end
  


end

if oe == 3 && type == 2
    quad = zeros(1,4);
    for i = 1:cp
    if i == 1
        l = 0;
    end
    [nodeCord  QuadNode, ln, last] = meshConc1(part(i,:), quad, node, elm, cp, count(i), count(i), lastnode, l);
    print(QuadNode,nodeCord);
    quad = QuadNode;
    node = nodeCord;
    lastnode = ln;
    nodeCord = node;
    l = last
    [rq cq] = size(QuadNode);
    [rn cn] = size(node);
    fprintf('ne = %d, total = %d', ne, total);
    if i == cp && ne == total
    [nodecord  QuadNode] = mesh3S(quad, node);
    print(QuadNode,nodecord);
    drawDomain(QuadNode,nodecord);
    [rq cq] = size(QuadNode);
    [rn cn] = size(node);
    end
    
    if i == cp && ne ~ total
        
    [nodeCord  QuadNode] = meshelm(quad, node, ne, n);
    print(QuadNode,nodeCord);
    quad = QuadNode;
    node = nodeCord; 
    lastnode = ln;
    [nodecord  QuadNode] = mesh3S(quad, node);
    print(QuadNode,nodecord);
    drawDomain(QuadNode, nodecord);
    [rq cq] = size(QuadNode);
    [rn cn] = size(nodecord);
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
                        %---Meshing Process (For Convex Polygon): For Case 1------%
                        
 function [node quadTri] = meshConv1(quad, node, ne, n)

if ne == 1
    return;
    
else
    quadTri = zeros(1,3)
    x = [0 0]
    for i = 1:n
        x = x+node(i,:);
    end
    x = x./n;
    node(i+1, :) = x;
   CentNode = i+1
   LastNode = CentNode +1
    for j = 1:n
        m = j+1;
        if j==n
            m=1;
        end
        node(LastNode, :) = (node(j,:) + node(m,:))/2;
        mid(j) = LastNode;
       LastNode = LastNode + 1;
    end
    for j = 1:n
        m = j+1;
        if j==n
            m=1;
        end
         quad(j,:) = [mid(j) m  mid(m) CentNode]
    end
     quad
    [noQuad , qr] = size(quad);
    [LastNode qn] = size(node); 
    elementsFlag = 1
    for i = 1:noQuad  
        x = [0 0] 
        for j = 1:qr
            x = (x+node(quad(i,j),:))
        end
        x = x/qr
        node(LastNode+1,:) = x
        CentNode(i) = LastNode + 1 
        LastNode = LastNode + 1
   for j = 1:4
        m = j+1
        if j == 4
            m = 1
        end
        elements(elementsFlag,:)  = [quad(i,j) quad(i,m) CentNode(i)]
        elementsFlag = elementsFlag + 1;
    end
    end  
    quadTri = elements
end

if ne == 4*n
    return;
    
else
 for l=1:10
        [noTri g] = size(quadTri);
        [lastNode g] = size(node);
       
        if noTri == ne
            break;
        end
       
       % tempNode = zeros(1,4);
        for i=1:noTri
             c = [0 0];
            for j=1:3
                m = j+1;
                if j==3
                    m=1;
                end
                
 x = (node(quadTri(i,j),:)+node(quadTri(i,m),:))/2;


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
 
    elements(elementsFlag, :) = [quadTri(i,1) tempNode(1) tempNode(3)];
    elements(elementsFlag+1, :) = [tempNode(1) quadTri(i,2) tempNode(2)];
    elements(elementsFlag+2, :) = [tempNode(3) tempNode(2) quadTri(i,3)]; 
    elements(elementsFlag+3, :) = [tempNode(1) tempNode(2) tempNode(3)];
    elementsFlag = elementsFlag + 4;
            
        end    
 quadTri = elements;

 end
end
end                       



                       %---Meshing Process (For Concave Polygon): For Case 1------% 
                       
function [node quad ln last] = meshConc1(part, quad, node, e, cp, ne, n, lastnode, last)

if e == cp
    quad(last+1,:) = part;
    last = last+1;
    ln = lastnode;
    node = node;
    return
end
    x = [0 0]
    fprintf ('ne = %d, n = %d last = %d', ne, n, last);
    for i = 1:n
        x = x+node(part(1,i),:);
    end
    x = x./n;
    node(lastnode, :) = x;
   CentNode = lastnode;
   LastNode = CentNode +1;
    for j = 1:n
        m = j+1;
        if j==n
            m=1;
        end
       x = (node(part(1,j),:) + node(part(1,m),:))/2;
        tmp = find(ismember(node,x,'rows'));
        if tmp<=lastnode
            mid(j) = tmp;
        else
        node(LastNode, :) = x;  
        mid(j) = LastNode;
       LastNode = LastNode + 1;
        end
    end
    for j = 1:n
        m = j+1;
        if j==n
            m=1;
        end
         quadtemp(j,:) = [mid(j) part(1,m)  mid(m) CentNode];
         ln = LastNode;
    end
    temp = 1;
    for time = last+1:last+ne
        quad(time,:) = quadtemp(temp,:)
        temp = temp+1;
        
    end
 last = time;

if ne == n
    return;
end
end

         %---Meshing Process: meshing each subdomain created from the problem domain ------%
         
function [node  quad] =  meshelm(quad, node, ne, n)
  for l=1:10
      fprintf ('ne in meshelm func = %d', ne);
        [noQuad g] = size(quad);
        [lastNode h] = size(node);
       
        if noQuad == ne
            break;
        end
        elementsFlag = 1;
        tempNode = zeros(1,4);
        for i=1:noQuad
             c = [0 0];
            for j=1:4
                m = j+1;
                if j==4
                    m=1;
                end
                
 x(1) = (node(quad(i,j),1)+node(quad(i,m),1))/2;
 
 x(2) = (node(quad(i,j),2)+node(quad(i,m),2))/2;
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
 
   c = c./4;
   temp = find(ismember(node,c,'rows'));
    m = j+1;
if temp<=lastNode
    tempNode(m) = temp;
else
    node(lastNode+1, :) = c;
    lastNode = lastNode + 1;
    tempNode(m) = lastNode;
end


    elements(elementsFlag, :) = [quad(i,1) tempNode(1) tempNode(5) tempNode(4)];
    elements(elementsFlag+1, :) = [tempNode(1) quad(i,2) tempNode(2) tempNode(5)];
    elements(elementsFlag+2, :) = [ tempNode(2) quad(i,3) tempNode(3) tempNode(5)];
    elements(elementsFlag+3, :) = [tempNode(3) quad(i, 4) tempNode(4) tempNode(5)];
    elementsFlag = elementsFlag + 4;
            
        end    
 quad = elements;

    end

end

%---Meshing Process: For Case 2------% 

function [node quad] = mesh2L(quad,node)
[ne tem] = size(quad);
[l temp2] = size(node);
for e = 1:ne
    temp = (node(quad(e,1), :) + node(quad(e,2),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
    temp = (node(quad(e,2), :) + node(quad(e,3),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end

     temp = (node(quad(e,3), :) + node(quad(e,4),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end

     temp = (node(quad(e,4), :) + node(quad(e,1),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      node(l+1,:) = (node(quad(e,1), :) + node(quad(e,2),:) + node(quad(e,3),:)+ node(quad(e,4),:))/4;
      l = l+1;
      q5 = l;
    
     
    
    newquad(e,:) = [quad(e,1) q1 quad(e,2) q2 quad(e,3) q3 quad(e,4) q4 q5];
    
end
quad = newquad;
return;

end


function [node quad] = mesh2S(quad,node)
[ne tem] = size(quad)
[l temp2] = size(node)
for e = 1:ne
    temp = (node(quad(e,1), :) + node(quad(e,2),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
    temp = (node(quad(e,2), :) + node(quad(e,3),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end

     temp = (node(quad(e,3), :) + node(quad(e,4),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end

     temp = (node(quad(e,4), :) + node(quad(e,1),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
    newquad(e,:) = [quad(e,1) q1 quad(e,2) q2 quad(e,3) q3 quad(e,4) q4];
end
quad = newquad;
return;

end

                        %---Meshing Process: For Case 3------% 
                        
                        
function [node quad] = mesh3L(quad,node)
[ne tem] = size(quad);
[l temp2] = size(node);
for e = 1:ne
    temp = (2*node(quad(e,1), :) + node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
     temp = (node(quad(e,1), :) + 2*node(quad(e,2),:))/3;
     q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
    temp = (2 * node(quad(e,2), :) + node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
 
    temp = (node(quad(e,2), :) + 2*node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
     temp = (2*node(quad(e,3), :) + node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
      temp = (node(quad(e,3), :) + 2* node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end

     temp = (2*node(quad(e,4), :) + node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
     temp = (node(quad(e,4), :) + 2*node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
 
      node(l+1,:) = (2*node(q1, :) + node(q6,:))/3;
      l = l+1;
      q9 = l;
      
      node(l+1,:) = (2*node(q2, :) + node(q5,:))/3;
      l = l+1;
      q10 = l;
     
      node(l+1,:) = (node(q2, :) + 2*node(q5,:))/3;
      l = l+1;
      q11 = l;
      
      node(l+1,:) = (node(q1, :) + 2*node(q6,:))/3;
      l = l+1;
      q12 = l;
    
    
    newquad(e,:) = [quad(e,1) q1 q2 quad(e,2) q3 q4 quad(e,3) q5 q6 quad(e,4) q7 q8 q9 q10 q11 q12];
end
quad = newquad;
return;

end
                        
                        
                        
function [node quad] = mesh3S(quad,node)
[ne tem] = size(quad);
[l temp2] = size(node);
for e = 1:ne
    temp = (2*node(quad(e,1), :) + node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
     temp = (node(quad(e,1), :) + 2*node(quad(e,2),:))/3;
     q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
    temp = (2 * node(quad(e,2), :) + node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
 
    temp = (node(quad(e,2), :) + 2*node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
     temp = (2*node(quad(e,3), :) + node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
      temp = (node(quad(e,3), :) + 2* node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end

     temp = (2*node(quad(e,4), :) + node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    
     temp = (node(quad(e,4), :) + 2*node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
    newquad(e,:) = [quad(e,1) q1 q2 quad(e,2) q3 q4 quad(e,3) q5 q6 quad(e,4) q7 q8];
end
quad = newquad;
return;

end

  %---Output 1: Prints the meshed elements and nodal coordinates of each nodes within these elements ------%


function print(quad, node)
[m n] = size(quad);
[k l] = size(node);
disp('Global Nodes (element wise): ');

for i=1:m
    fprintf('\n[%2i]', i);
    for j = 1:n
        fprintf('%3i ',quad(i,j));
    end
end

fprintf('\nNodal Co-ordinates: \n\n');
fprintf('Node NO.:  abscissa  ordinate\n');

for i=1:k
    fprintf('  %3i : %0.3f    %0.3f\n', i, node(i,1), node(i,2));

end
end


    %---Output 2: Plots the figure of the meshed domain------%

     function drawDomain(q, c)
        [k l] = size(q);
        [r s] = size(c);
        
 
              
        xmin = min(c(:,1)) ;
        ymin = min(c(:,2));
        xmax = max(c(:,1)) ;
        ymax = max(c(:,2));
        
        figure
        xlim([xmin-0.5, xmax+0.5]); ylim([ymin-0.5, ymax+0.5]);  
        hold on ; 
        
        for i=1:k
            cx = sum(c(q(i,:), 1))/l
            cy = sum(c(q(i,:), 2))/l
            if l<16
            for j = 1:l
                m = j+1
                if l == 9
                x1 = c(q(i,m),1)
                y1 = c(q(i,m),2)
                
                value1 = abs(cx - x1)
                value2 = abs(cy - y1)
                if abs([cx cy] - [x1 y1])<=10^-10
                    m = 1
                x = [c(q(i,j),1) c(q(i,m),1)];
                y = [c(q(i,j),2) c(q(i,m),2)];
                plot(x,y);
                hold on; 
                break; 
                end
                end
           
                if j==l && m ~ 1
                    m=1;
                end
                
                x = [c(q(i,j),1) c(q(i,m),1)];
                y = [c(q(i,j),2) c(q(i,m),2)];
               
                plot(x,y);
                hold on 
            end
            else
            for j = 1:12
                
                      m = j+1;
             
                      if j == 12
                          m =1;
                      end
                
                x = [c(q(i,j),1) c(q(i,m),1)];
                y = [c(q(i,j),2) c(q(i,m),2)];
               
                plot(x,y);
                hold on
            end
             end   
       
           
          % text(cx-0.06,cy,'(');
           %  hold on
           % text(cx-0.04,cy,int2str(i));
           % hold on
           % text(cx-0.02,cy,')');
           % hold on
            
        end
        
       for i = 1:k
            for  j = 1:l
                if l == 4 || l == 8 || l == 12
                           plot(c(q(i,j),1),c(q(i,j),2), 'Marker', 'O', 'MarkerFaceColor', 'red');
                           hold on
                           text(c(q(i,j),1)-0.05, c(q(i,j),2)-0.05, int2str(q(i,j)));
                           hold on
                elseif l == 9
                    if j == 9
                         plot(c(q(i,j),1),c(q(i,j),2), 'Marker', 'O', 'MarkerFaceColor', 'blue');
                         hold on
                         text(c(q(i,j),1)-0.05, c(q(i,j),2)-0.05, int2str(q(i,j)));
                          hold on
                    else
                        plot(c(q(i,j),1),c(q(i,j),2), 'Marker', 'O', 'MarkerFaceColor', 'red');
                        hold on
                        text(c(q(i,j),1)-0.05, c(q(i,j),2)-0.05, int2str(q(i,j)));
                        hold on
                    end
                
                elseif l == 16
                        if j>=13
                             plot(c(q(i,j),1),c(q(i,j),2), 'Marker', 'O', 'MarkerFaceColor', 'blue');
                             hold on
                                 text(c(q(i,j),1)-0.05, c(q(i,j),2)-0.05, int2str(q(i,j)));
                               hold on
                        else
                              plot(c(q(i,j),1),c(q(i,j),2), 'Marker', 'O', 'MarkerFaceColor', 'red');
                              hold on
                              text(c(q(i,j),1)-0.05, c(q(i,j),2)-0.05, int2str(q(i,j)));
                              hold on
                        end  
         
            
                end
            end
        end
        xlabel('x axis');
        ylabel('y axis');

        
     end