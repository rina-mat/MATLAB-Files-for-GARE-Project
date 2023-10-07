function Assembled_Program_for_curvature_sided_polygonal_domain_trngl
syms x
format short
n =input('Number of vertices? ');

fprintf('Instruction: 1) Enter the nodal co-ordinates in an ascending order/descending order \n \n');
fprintf('e.g. [-1 0] -> [0 0] -> [1 0] or in reverse \n \t \t \n');
 
  %%% INPUT : Domain Co-ordinates, No. of elements of required order, Type
  %%% of element Lagrange/Serendipity 

fprintf('Enter the co-ordinates: \n');

for i = 1:n
    fprintf('Node %i: ', i);
    node(i, :) = input(''); 
    poly(1,i) = i;
end

if n == 3
   tri = poly;
else 
    tri = zeros(1,3);
    
end

oe = input('Enter order of element(1 or 2 or 3):');

type_of_domain = input('Does the domain have curve side/sides? Yes: Not all sides curved = 1, All sides Curved = 2, No = 0;    ');
if type_of_domain ~= 1 && type_of_domain ~= 0 && type_of_domain ~= 2
    fprintf('WARNING!! Please enter a valid input(0/1/2);');
    type_of_domain = input(' ');
end
if type_of_domain == 1 
curve_type = input('Type of curvatures:  \n  curve of same function = 1, curve of different function = 2;       ');
if curve_type ~= 1 && curve_type ~= 2
    fprintf('WARNING!! Please press 1 or 2;');
    curve_type = input(' '); 
end
number_of_curve_sides = input('How many curved sides are there? ');
ref = input('Refinement No. :  ');
ne = 4^(ref-1)*n;

if curve_type == 1
  
    if number_of_curve_sides == 1
            gg = input('Enter function of the curvature:  '); 
            Tri_Meshing_One_Curve_Side(tri,node,ne,n,gg,oe); 
            return;
    else
              gg = input('Enter function of the curvature:  '); 
            Tri_Meshing_more_than_one_Curve_Sides(tri,node,ne,n,gg,oe);
            return;
    end
    
else
    Num_of_fn = number_of_curve_sides;
for i = 1:Num_of_fn
fprintf('Expression of curve %d: ', i);
curve{i} = input('');
end
TriMeshing_With_Boundaries_Having_Different_Curve_Fn(tri,node,ne,n,curve,oe);
return;
end

elseif type_of_domain == 2
    ref = input('Refinement No. :  ');
    ne = 4^(ref-1)*n;
    Num_of_fn = n;
    for i = 1:Num_of_fn
fprintf('Expression of curve %d: ', i);
curve{i} = input('');

    end
TriMeshing_boundary_enclosed_by_curves(tri,node,ne,n,curve,oe)
else
    
    %initial domain(before meshing)
        xmin = min(node(:,1)) ;
        ymin = min(node(:,2));
        xmax = max(node(:,1)) ;
        ymax = max(node(:,2));
        
        figure
        xlim([xmin-0.5, xmax+0.5]); ylim([ymin-0.5, ymax+0.5]);  
        hold on ; 
        for i = 1:n
            j = i+1;
            if i == n
                j =1;
            end
            x = [node(i,1) node(j,1)];
            y = [node(i,2) node(j,2)];
            plot(x,y);
            hold on;
            plot(node(i,1), node(i,2), 'o', 'MarkerSize', 4,'MarkerFaceColor', 'red');
        end
        title('The domain to be meshed');
        xlabel('x axis');
        ylabel('y axis');
        
       
  Tri_Meshing_Linear_sided_Convex_Polygon(tri,node,n,oe,type);
  return;
end
 
end 

    
    
  
function Tri_Meshing_Linear_sided_Convex_Polygon(tri,node,n,oe)

maxref = input('Enter Number of maximum refinements (Put any integer value (1 to 4 is preferable)): ');

for i = 1:maxref
    ne(i) = 4^(i-1)*n;
end

for i = 1: maxref 
    if n == 3
    tri = poly;
else 
    tri = zeros(1,4);
end
if oe == 1
    
 
    [nodeCord, TriNode] = mesh1_lin(tri, node, ne(i), n);
    fprintf('For Refinement No.: %i', i);
    print(TriNode,nodeCord);
    tri = TriNode;
    node = nodeCord;
    drawDomain_linear(tri, node);
  
end

%%% Quadratic element of Lagrange Type %%%
  if oe == 2 

    [nodeCord,  TriNode] = mesh1_lin(tri, node, ne(i), n); 
    tri = TriNode;
    node = nodeCord;
   [nodecord,  TriNode] = mesh2_lin(tri, node);
    fprintf('For Refinement No.: %i', i);
    print(TriNode,nodecord);
    drawDomain_linear(TriNode,nodecord);
  end


if oe == 3 
    [nodeCord,  TriNode] = mesh1_lin(tri, node, ne(i), n); 
    tri = TriNode;
    node = nodeCord;
   [nodecord,  TriNode] = mesh3_lin(tri, node);
    fprintf('For Refinement No.: %i', i);
    print(TriNode,nodecord);
   drawDomain_linear(TriNode,nodecord);
  
end

end


end

function [node, tri] = mesh1_lin(tri, node, ne, n)

if ne == 1
    return;
    
else
    x = [0 0]
    for i = 1:n
        x = x+node(i,:);
    end
    x = x./n;
    node(i+1, :) = x;
   CentNode = i+1
   lastNode = CentNode 
   
    for j = 1:n
        m = j+1;
        if j==n
            m=1;
        end
         tri(j,:) = [j m CentNode];
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

function [node, quad] = mesh2_lin(tri,node) 
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

function [node, tri] = mesh3_lin(tri,node)
[ne tem] = size(tri);
[l temp2] = size(node);
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

%%% END OF PROGRAM FOR MESHING LINEAR CONVEX POLYGONAL DOMAIN %%%

 function Tri_Meshing_One_Curve_Side(tri,node,ne,n,gg,oe)
 
 cv = input('Enter boundary nodes which are in the curvature part: ');
icv = cv; %restoring initial array
[cvr, cvc] = size(cv)
count = cvc+1;
    %%%initial domain (before meshing)
        
             figure  
 
            cmax = max(node(icv(1),1), node(icv(cvc), 1));
            cmin = min(node(icv(1),1), node(icv(cvc), 1));
            fplot(gg,[cmin cmax], 'b');
            hold on;
       
        
for j = 1:n
      for i = 1:count-1
            if cv(i) == j 
                position(j) = i
            end
      end
end
 for j = 1:n
     plot(node(j,1), node(j,2), 'o', 'MarkerSize', 4,'MarkerFaceColor', 'red');
     hold on;
        m = j+1;
        if j==n
            m=1;
        end
       
        if (any(cv(:) == j) && any(cv(:) == m)) && (position(m)>position(j))
            continue;
        else
            x = [node(j,1) node(m,1)];
            y = [node(j,2) node(m,2)];
            plot(x,y);
            hold on;
        end
         
 end
   title ('The domain to be meshed')
   xlabel('x axis');
    ylabel('y axis');
 
 if oe == 1
    
 
    [nodeCord  TriNode cv count] = mesh_os1(tri, node,  gg,  ne, n, cv, count);

    print(TriNode,nodeCord);
    tri = TriNode;
    node = nodeCord;
    nodeCord = node;
    drawDomain_os(tri, node, icv, gg);

end

%%% Quadratic element of Lagrange Type %%%
  if oe == 2

    [nodeCord  TriNode cv count] = mesh_os1(tri, node, gg,  ne, n, cv, count); 
    tri = TriNode;
    node = nodeCord;
   [nodecord  TriNode cv] = mesh_os2(tri, node, gg, cv, count);
    print(TriNode,nodecord);
    drawDomain_os(TriNode,nodecord, icv, gg);
  end



%%% Cubic element of Lagrange Type %%%
if oe == 3
     [nodeCord  TriNode cv count] = mesh_os1(tri, node, gg,  ne, n, cv, count); 
    tri = TriNode;
    node = nodeCord;
 
   [nodecord  TriNode cv] = mesh_os3(tri, node, gg, cv, count);
  
    print(TriNode,nodecord);
   drawDomain_os(TriNode,nodecord, icv, gg);   
end
 end

 function [node tri cv count] = mesh_os1(tri, node, gg, ne, n, cv, count)

lastNode = n+1;

for j = 1:n
      for i = 1:count-1
            if cv(i) == j 
                pos(j) = i
            end
      end
end
    
     x = [0 0];
    for i = 1:lastNode-1
        x = x+node(i,:);
    end
    x = x/(lastNode-1);
    node(lastNode, :) = x;
    CentNode = lastNode;
    lastNode = lastNode+1;
  
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
    
        for i=1:noTri
            for j=1:3
                m = j+1;
                if j==3
                    m=1;
                end
  if any(cv(:) == tri(i,j)) && any(cv(:) == tri(i,m))
   fprintf('Functional midvalue of  [%0.3f %0.3f] and [%0.3f %0.3f]', node(tri(i,j),1), node(tri(i,j),2), node(tri(i,m),1), node(tri(i,m),2));
   x(1) = (node(tri(i,j),1)+node(tri(i,m),1))/2;
 
   x(2) = gg(x(1));
   x = [x(1) x(2)];
   
   temp = find(ismember(node,x,'rows'));
   
if temp<=lastNode
   tempNode(j) = temp;
    
   else
    node(lastNode+1, :) = x; 
    lastNode = lastNode + 1;
    cv(count) = lastNode; 
    count = count+1;
    tempNode(j) = lastNode;
   
   
end
  else
     x(1) = (node(tri(i,j),1)+node(tri(i,m),1))/2;
 
     x(2) = (node(tri(i,j),2)+node(tri(i,m),2))/2;
     
     x = [x(1) x(2)];
     temp = find(ismember(node,x,'rows'));



if temp<=lastNode
    tempNode(j) = temp;
     
else
    node(lastNode + 1, :) = x;
    lastNode = lastNode + 1;   
    tempNode(j) = lastNode;
   
end

            end
 end
    elements(elementsFlag, :) = [tri(i,1) tempNode(1) tempNode(3)];
   elements(elementsFlag+1, :) = [tempNode(1) tri(i,2) tempNode(2)];
   elements(elementsFlag+2, :) = [tempNode(3) tempNode(2) tri(i,3)]; 
   elements(elementsFlag+3, :) = [tempNode(1) tempNode(2) tempNode(3)];
    elementsFlag = elementsFlag + 4;
 end       
 tri = elements;
 cv
        
end
end
end



function [node quad cv] = mesh_os2(quad,node, gg, cv, count)
[ne tem] = size(quad);
[l temp2] = size(node);
for e = 1:ne
      if any(cv(:) == quad(e,1)) && any(cv(:) == quad(e,2))
           x(1) = (node(quad(e,1),1)+node(quad(e,2),1))/2;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(count) = l;
        count  = count+1;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + node(quad(e,2),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
    
      
      if any(cv(:) == quad(e,2)) && any(cv(:) == quad(e,3))
           x(1) = (node(quad(e,2),1)+node(quad(e,3),1))/2;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(count) = l;
        count  = count+1;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + node(quad(e,3),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
      
      if any(cv(:) == quad(e,3)) && any(cv(:) == quad(e,4))
           x(1) = (node(quad(e,3),1)+node(quad(e,4),1))/2;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(count) = l;
        count  = count+1;
        node(l,:) = temp;
    end
      else
     temp = (node(quad(e,3), :) + node(quad(e,4),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
      
      if any(cv(:) == quad(e,4)) && any(cv(:) == quad(e,1))
           x(1) = (node(quad(e,4),1)+node(quad(e,1),1))/2;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(count) = l;
        count  = count+1;
        node(l,:) = temp;
    end
      else
     
     temp = (node(quad(e,4), :) + node(quad(e,1),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
    
    node(l+1,:) = (node(quad(e,1), :) + node(quad(e,2),:) + node(quad(e,3),:)+ node(quad(e,4),:))/4;
      l = l+1;
      q5 = l;
    
     
    
    newquad(e,:) = [quad(e,1) q1 quad(e,2) q2 quad(e,3) q3 quad(e,4) q4 q5];
end
quad = newquad;
return;
end




function [node quad cv] = mesh_os3(quad,node, gg, cv, count)
[ne tem] = size(quad);
[l temp2] = size(node);
for e = 1:ne
    
         if any(cv(:) == quad(e,1)) && any(cv(:) == quad(e,2))
           x(1) = (2*node(quad(e,1), 1) + node(quad(e,2),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(count) = l;
        count  = count+1;
        node(l,:) = temp;
    end
         else
    temp = (2*node(quad(e,1), :) + node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end
    
    if any(cv(:) == quad(e,1)) && any(cv(:) == quad(e,2))
           x(1) = (node(quad(e,1), 1) + 2*node(quad(e,2),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(count) = l;
        count  = count+1;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + 2*node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end
    
            if any(cv(:) == quad(e,2)) && any(cv(:) == quad(e,3))
           x(1) = (2*node(quad(e,2), 1) + node(quad(e,3),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(count) = l;
        count  = count+1;
        node(l,:) = temp;
    end
     else
    temp = (2*node(quad(e,2), :) + node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
            end
            
      if any(cv(:) == quad(e,2)) && any(cv(:) == quad(e,3))
          x(1) = (node(quad(e,2), 1) + 2*node(quad(e,3),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(count) = l;
        count  = count+1;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + 2*node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
            end      
   
        if any(cv(:) == quad(e,3)) && any(cv(:) == quad(e,4))
           x(1) = (2*node(quad(e,3), 1) + node(quad(e,4),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        cv(count) = l;
        count  = count+1;
        node(l,:) = temp;
    end
        else
    temp = (2*node(quad(e,3), :) + node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
        end
          if any(cv(:) == quad(e,3)) && any(cv(:) == quad(e,4))
          x(1) = (node(quad(e,3), 1) + 2*node(quad(e,4),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        cv(count) = l;
        count  = count+1;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,3), :) + 2*node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end 
    
        if any(cv(:) == quad(e,4)) && any(cv(:) == quad(e,1))
           x(1) = (2*node(quad(e,4), 1) + node(quad(e,1),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        cv(count) = l;
        count  = count+1;
        node(l,:) = temp;
    end
        else
    temp = (2*node(quad(e,4), :) + node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
        end
        if any(cv(:) == quad(e,4)) && any(cv(:) == quad(e,1))
          x(1) = (node(quad(e,4), 1) + 2*node(quad(e,1),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        cv(count) = l;
        count  = count+1;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,4), :) + 2*node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
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


function Tri_Meshing_more_than_one_Curve_Sides(tri,node,ne,n,gg,oe)

cv = input('Enter boundary nodes which are in the curvature part: ');
[cvr cvc] = size(cv)
if oe == 1
    
 
    [nodeCord,  TriNode,  cv] = mesh1(tri, node,  gg,  ne, n, cv);

    print(TriNode,nodeCord);
    tri = TriNode;
    node = nodeCord;
   drawDomain(tri, node, cv, gg);
end

if oe == 2  

    [nodeCord,  TriNode, cv] = mesh1(tri, node, gg,  ne, n, cv); 
    tri = TriNode;
    node = nodeCord;
   [nodecord,  TriNode] = mesh2(tri, node, gg, cv);
    print(TriNode,nodecord);
    drawDomain(TriNode,nodecord, cv, gg);
  end

if oe == 3
    [nodeCord,  TriNode, cv] = mesh1(tri, node, gg,  ne, n, cv); 
    tri = TriNode;
    node = nodeCord;
   [nodecord,  TriNode, cv] = mesh3(tri, node, gg, cv)
    print(TriNode,nodecord);
   drawDomain(TriNode,nodecord, cv, gg);   
end

end
 



function [node, tri, cv] = mesh1(tri, node, gg, ne, n, cv)
LastNode = n+1;

     x = [0 0];
    for i = 1:LastNode-1
        x = x+node(i,:);
    end
    x = x/(LastNode-1);
    node(LastNode, :) = x;
    CentNode = LastNode;
    LastNode = LastNode+1;
  
    for j = 1:n
        m = j+1;
        if j==n
            m=1;
        end
         tri(j,:) = [j m CentNode];
       
    end
    


if ne == n
    return;
    
else
    for l=1:10
        [noTri g] = size(tri)
        [lastNode g] = size(node)
       
        if noTri == ne
            break;
        end
        elementsFlag = 1;
        %tempNode = zeros(1,4);
        for i=1:noTri
             c = [0 0];
             [cvr cvc] = size(cv)
            for j=1:3
                m = j+1;
                if j==3
                    m=1;
                end
         [cj rj] = find(cv == tri(i,j))
         if isempty(find(cv == tri(i,j)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == tri(i,m))
         if isempty(find(cv == tri(i,m)))
             cm = -1; rm = -1;
         end
  if cj == cm
   x(1) = (node(tri(i,j),1)+node(tri(i,m),1))/2;
 
   x(2) = gg(x(1));
   x = [x(1) x(2)];
   
   temp = find(ismember(node,x,'rows'));
   
if temp<=lastNode
   tempNode(j) = temp;
    c = c + node(temp, :);
    
   else
    node(lastNode+1, :) = x;
    lastNode = lastNode + 1;
    cv(cj, cvc+1) = lastNode;
    tempNode(j) = lastNode;
    c = c + node(lastNode, :);
end
  else
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
        xmin = min(node(:,1)) ;
        ymin = min(node(:,2));
        xmax = max(node(:,1)) ;
        ymax = max(node(:,2));
        
        figure
        xlim([xmin-0.5, xmax+0.5]); ylim([ymin-0.5, ymax+0.5]);  
        hold on ; 
       
        for i=1:cvr
            cmax = max(node(cv(i,1),1), node(cv(i,2), 1));
            cmin = min(node(cv(i,1),1), node(cv(i,2), 1));
            fplot(gg,[cmin cmax], 'b');
        end
for j = 1:n
        m = j+1;
        if j==n
            m=1;
        end
         [cj, rj] = find(cv == j)
         if isempty(find(cv == j))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == m)
         if isempty(find(cv == m))
             cm = -1; rm = -1;
         end
         
         if cj == cm
             continue
         else
             x = [node(j,1) node(m,1)];
             y = [node(j,2) node(m,2)];
             plot(x,y);
             hold on;
         end
          plot(node(j,1), node(j,2), 'o', 'MarkerSize', 4,'MarkerFaceColor', 'red');
 end
         title ('The domain to be meshed')
         xlabel('x axis');
         ylabel('y axis');
end

function [node, quad, cv] = mesh2L(quad,node, gg, cv)
[ne, tem] = size(quad);
[l, temp2] = size(node);
[cvr, cvc] = size(cv)
for e = 1:ne
        [cj, rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,1),1)+node(quad(e,2),1))/2;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + node(quad(e,2),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
    
      
         [cj, rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,2),1)+node(quad(e,3),1))/2;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + node(quad(e,3),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
      
      [cj, rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,3),1)+node(quad(e,4),1))/2;
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
     temp = (node(quad(e,3), :) + node(quad(e,4),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
     end
      
        [cj, rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,4),1)+node(quad(e,1),1))/2;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
     
     temp = (node(quad(e,4), :) + node(quad(e,1),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
    
    node(l+1,:) = (node(quad(e,1), :) + node(quad(e,2),:) + node(quad(e,3),:)+ node(quad(e,4),:))/4;
      l = l+1;
      q5 = l;
    
     
    
    newquad(e,:) = [quad(e,1) q1 quad(e,2) q2 quad(e,3) q3 quad(e,4) q4 q5];
end
quad = newquad;
return;
end


function [node, quad, cv] = mesh2S(quad,node, gg, cv)
[ne, tem] = size(quad);
[l, temp2] = size(node);
[cvr, cvc] = size(cv)
for e = 1:ne
        [cj, rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,1),1)+node(quad(e,2),1))/2;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + node(quad(e,2),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
    
      
         [cj, rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,2),1)+node(quad(e,3),1))/2;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + node(quad(e,3),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
      
      [cj, rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,3),1)+node(quad(e,4),1))/2;
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
     temp = (node(quad(e,3), :) + node(quad(e,4),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
     end
      
        [cj, rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,4),1)+node(quad(e,1),1))/2;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
     
     temp = (node(quad(e,4), :) + node(quad(e,1),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
     end
     
    
    newquad(e,:) = [quad(e,1) q1 quad(e,2) q2 quad(e,3) q3 quad(e,4) q4];
end
quad = newquad;
return;

end


function [node, quad, cv] = mesh3L(quad,node, gg, cv)
[ne, tem] = size(quad);
[l, temp2] = size(node);
[cvr, cvc] = size(cv);
for e = 1:ne
    
           [cj, rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (2*node(quad(e,1), 1) + node(quad(e,2),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
         else
    temp = (2*node(quad(e,1), :) + node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end
    
  [cj, rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,1), 1) + 2*node(quad(e,2),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + 2*node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end
    
          [cj, rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (2*node(quad(e,2), 1) + node(quad(e,3),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
     else
    temp = (2*node(quad(e,2), :) + node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
            end
            
       [cj, rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if cj == cm
          x(1) = (node(quad(e,2), 1) + 2*node(quad(e,3),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + 2*node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
            end      
    [cj, rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (2*node(quad(e,3), 1) + node(quad(e,4),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
        else
    temp = (2*node(quad(e,3), :) + node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
        end
      [cj, rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if cj == cm
          x(1) = (node(quad(e,3), 1) + 2*node(quad(e,4),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,3), :) + 2*node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end 
    
      [cj, rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (2*node(quad(e,4), 1) + node(quad(e,1),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
        else
    temp = (2*node(quad(e,4), :) + node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
        end
  [cj, rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if cj == cm
          x(1) = (node(quad(e,4), 1) + 2*node(quad(e,1),1))/3;
 
           x(2) = gg(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,4), :) + 2*node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
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


function TriMeshing_With_Boundaries_Having_Different_Curve_Fn(tri,node,ne,n,curve,oe)

cv = input('Enter boundary nodes which are in the curvature part: ');

if oe == 1
    
 
    [nodeCord  TriNode cv] = mesh_DCF_1(tri, node,  curve,  ne, n, cv);

    print(TriNode,nodeCord);
    tri = TriNode;
    node = nodeCord;
   drawDomain_DCF(tri, node, cv, curve);
    end

%%% Quadratic element of Lagrange Type %%%
if oe == 2

    [nodeCord  TriNode cv] = mesh_DCF_1(tri, node,  curve,  ne, n, cv); 
    tri = TriNode;
    node = nodeCord;
   [nodecord  TriNode cv] = mesh_DCF_2L(tri, node,  curve, cv);
    print(TriNode,nodecord);
    drawDomain_DCF(TriNode,nodecord, cv,  curve);
    [rq cq] = size(TriNode);
    [rn cn] = size(nodecord);
  end

%%% Cubic element of Lagrange Type %%%
if oe == 3
    [nodeCord  TriNode cv] = mesh_DCF_1(tri, node,  curve,  ne, n, cv); 
    tri = TriNode;
    node = nodeCord;
   [nodecord  TriNode cv] = mesh_DCF_3L(tri, node,  curve, cv)
    print(TriNode,nodecord);
   drawDomain_DCF(TriNode,nodecord, cv,  curve);
    [rq cq] = size(TriNode);
    [rn cn] = size(nodecord);
    
end
end

function [node tri cv] = mesh_DCF_1(tri, node, curve, ne, n, cv)
LastNode = n+1

 x = [0 0];
    for i = 1:LastNode-1
        x = x+node(i,:);
    end
    x = x/(LastNode-1);
    node(LastNode, :) = x;
    CentNode = LastNode;
    LastNode = LastNode+1;
  
    for j = 1:n
        m = j+1;
        if j==n
            m=1;
        end
         tri(j,:) = [j m CentNode];
       
    end
 
if ne == n
    return;
    
else
    for l=1:10
        [noTri g] = size(tri)
        [lastNode g] = size(node)
       
        if noTri == ne
            break;
        end
        elementsFlag = 1;
        %tempNode = zeros(1,4);
        for i=1:noTri
             c = [0 0];
             [cvr cvc] = size(cv)
            for j=1:3
                m = j+1;
                if j==3
                    m=1;
                end
         [cj rj] = find(cv == tri(i,j))
         if isempty(find(cv == tri(i,j)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == tri(i,m))
         if isempty(find(cv == tri(i,m)))
             cm = -1; rm = -1;
         end
  if cj == cm
   x(1) = (node(tri(i,j),1)+node(tri(i,m),1))/2;
 
   x(2) = curve{cj}(x(1));
   x = [x(1) x(2)];
   
   temp = find(ismember(node,x,'rows'));
   
if temp<=lastNode
   tempNode(j) = temp;
    
   else
    node(lastNode+1, :) = x;
    lastNode = lastNode + 1;
    cv(cj, cvc+1) = lastNode;
    tempNode(j) = lastNode;
end
  else
     x(1) = (node(tri(i,j),1)+node(tri(i,m),1))/2;
     x(2) = (node(tri(i,j),2)+node(tri(i,m),2))/2;
     x = [x(1) x(2)];
     temp = find(ismember(node,x,'rows'));

if temp<=lastNode
    tempNode(j) = temp;    
else
    node(lastNode+1, :) = x;
    lastNode = lastNode + 1;
    tempNode(j) = lastNode;
  
end
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
%%%initial domain (before meshing)
        xmin = min(node(:,1)) ;
        ymin = min(node(:,2));
        xmax = max(node(:,1)) ;
        ymax = max(node(:,2));
        
        figure
        xlim([xmin-0.5, xmax+0.5]); ylim([ymin-0.5, ymax+0.5]);  
        hold on ; 
       
        for i=1:cvr
            cmax = max(node(cv(i,1),1), node(cv(i,2), 1));
            cmin = min(node(cv(i,1),1), node(cv(i,2), 1));
            fplot(curve{i},[cmin cmax], 'b');
        end
for j = 1:n
        plot(node(j,1), node(j,2), 'o', 'MarkerSize', 4,'MarkerFaceColor', 'red');
         hold on;
        m = j+1;
        if j==n
            m=1;
        end
         [cj, rj] = find(cv == j)
         if isempty(find(cv == j))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == m)
         if isempty(find(cv == m))
             cm = -1; rm = -1;
         end
         
         if cj == cm
             continue
         else
             x = [node(j,1) node(m,1)];
             y = [node(j,2) node(m,2)];
             plot(x,y);
             hold on;
         end
        
 end
         title ('The domain to be meshed')
         xlabel('x axis');
        ylabel('y axis');
end



function [node quad cv] = mesh_DCF_2L(quad,node,  curve, cv)
[ne tem] = size(quad);
[l temp2] = size(node);
[cvr cvc] = size(cv)
for e = 1:ne
        [cj rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,1),1)+node(quad(e,2),1))/2;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + node(quad(e,2),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
    
      
         [cj rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,2),1)+node(quad(e,3),1))/2;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + node(quad(e,3),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
      
      [cj rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,3),1)+node(quad(e,4),1))/2;
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
     temp = (node(quad(e,3), :) + node(quad(e,4),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
     end
      
        [cj rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,4),1)+node(quad(e,1),1))/2;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
     
     temp = (node(quad(e,4), :) + node(quad(e,1),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
    
    node(l+1,:) = (node(quad(e,1), :) + node(quad(e,2),:) + node(quad(e,3),:)+ node(quad(e,4),:))/4;
      l = l+1;
      q5 = l;
    
     
    
    newquad(e,:) = [quad(e,1) q1 quad(e,2) q2 quad(e,3) q3 quad(e,4) q4 q5];
end
quad = newquad;
return;
end


function [node quad cv] = mesh_DCF_2S(quad,node,  curve, cv)
[ne tem] = size(quad);
[l temp2] = size(node);
[cvr cvc] = size(cv)
for e = 1:ne
        [cj rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,1),1)+node(quad(e,2),1))/2;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + node(quad(e,2),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
    
      
         [cj rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,2),1)+node(quad(e,3),1))/2;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + node(quad(e,3),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
      
      [cj rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,3),1)+node(quad(e,4),1))/2;
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
     temp = (node(quad(e,3), :) + node(quad(e,4),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
     end
      
        [cj rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,4),1)+node(quad(e,1),1))/2;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
     
     temp = (node(quad(e,4), :) + node(quad(e,1),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
     end
     
    
    newquad(e,:) = [quad(e,1) q1 quad(e,2) q2 quad(e,3) q3 quad(e,4) q4];
end
quad = newquad;
return;

end


function [node quad cv] = mesh_DCF_3L(quad,node,  curve, cv)
cv
[ne tem] = size(quad);
[l temp2] = size(node);
[cvr cvc] = size(cv);
for e = 1:ne
    
        [cj, rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (2*node(quad(e,1), 1) + node(quad(e,2),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
         else
    temp = (2*node(quad(e,1), :) + node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end
    
  [cj, rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,1), 1) + 2*node(quad(e,2),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + 2*node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end
    
          [cj, rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (2*node(quad(e,2), 1) + node(quad(e,3),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
     else
    temp = (2*node(quad(e,2), :) + node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
            end
            
       [cj, rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if cj == cm
          x(1) = (node(quad(e,2), 1) + 2*node(quad(e,3),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + 2*node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
            end      
    [cj, rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (2*node(quad(e,3), 1) + node(quad(e,4),1))/3;
 
           x(2) = curve(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
        else
    temp = (2*node(quad(e,3), :) + node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
        end
      [cj, rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if cj == cm
          x(1) = (node(quad(e,3), 1) + 2*node(quad(e,4),1))/3;
 
           x(2) = curve(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,3), :) + 2*node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end 
    
      [cj, rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (2*node(quad(e,4), 1) + node(quad(e,1),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
        else
    temp = (2*node(quad(e,4), :) + node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
        end
  [cj, rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if cj == cm
          x(1) = (node(quad(e,4), 1) + 2*node(quad(e,1),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,4), :) + 2*node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
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


function [node quad cv] = mesh_DCF_3S(quad,node,  curve, cv)
[ne tem] = size(quad);
[l temp2] = size(node);
[cvr cvc] = size(cv);
for e = 1:ne
    
           [cj rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (2*node(quad(e,1), 1) + node(quad(e,2),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
         else
    temp = (2*node(quad(e,1), :) + node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end
    
  [cj rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2), 1))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (node(quad(e,1), 1) + 2*node(quad(e,2),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + 2*node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end
    
          [cj rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (2*node(quad(e,2), 1) + node(quad(e,3),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
     else
    temp = (2*node(quad(e,2), :) + node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
            end
            
       [cj rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if cj == cm
          x(1) = (node(quad(e,2), 1) + 2*node(quad(e,3),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + 2*node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
            end      
    [cj rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (2*node(quad(e,3), 1) + node(quad(e,4),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
        else
    temp = (2*node(quad(e,3), :) + node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
        end
      [cj rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if cj == cm
          x(1) = (node(quad(e,3), 1) + 2*node(quad(e,4),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,3), :) + 2*node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end 
    
      [cj rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if cj == cm
           x(1) = (2*node(quad(e,4), 1) + node(quad(e,1),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
        else
    temp = (2*node(quad(e,4), :) + node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
        end
  [cj rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if cj == cm
          x(1) = (node(quad(e,4), 1) + 2*node(quad(e,1),1))/3;
 
           x(2) = curve{cj}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        cv(cj,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,4), :) + 2*node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end 
    
    newquad(e,:) = [quad(e,1) q1 q2 quad(e,2) q3 q4 quad(e,3) q5 q6 quad(e,4) q7 q8];
    
end
quad = newquad;
return;
end
%%% Boundary Enclosed by Curves

function TriMeshing_boundary_enclosed_by_curves(tri,node,ne,n,curve,oe,type)
cv = input('End points of each curve boundary(minimum 3 points):   ');

if oe == 1
    
 
    [nodeCord  TriNode cv] = mesh_BEC_1(tri, node,  curve,  ne, n, cv);

    print(TriNode,nodeCord);
    tri = TriNode;
    node = nodeCord;
    [cvr cvc] = size(cv);
    %nodeCord = node;

   drawDomain_DCF(tri, node, cv, curve);
    [rq cq] = size(TriNode)
    [rn cn] = size(node)

end

%%% Quadratic element of Lagrange Type %%%
if oe == 2 && type == 1 

    [nodeCord  TriNode cv] = mesh_BEC_1(tri, node,  curve,  ne, n, cv); 
    tri = TriNode;
    node = nodeCord;
   [nodecord  TriNode cv] = mesh_BEC_2L(tri, node,  curve, cv);
    print(TriNode,nodecord);
    drawDomain_DCF(TriNode,nodecord, cv,  curve);
    [rq cq] = size(TriNode);
    [rn cn] = size(nodecord);
end

%%% Cubic element of Lagrange Type %%%
if oe == 3 && type == 1
    [nodeCord  TriNode cv] = mesh_BEC_1(tri, node,  curve,  ne, n, cv); 
    tri = TriNode;
    node = nodeCord;
   [nodecord  TriNode cv] = mesh_BEC_3L(tri, node,  curve, cv)
    print(TriNode,nodecord);
   drawDomain_DCF(TriNode,nodecord, cv,  curve);
    [rq cq] = size(TriNode);
    [rn cn] = size(nodecord);
    
end
end
function [node tri cv] = mesh_BEC_1(tri, node, curve, ne, n, cv)

LastNode = n+1;
[cvr cvc] = size(cv);

     x = [0 0];
    for i = 1:LastNode-1
        x = x+node(i,:);
    end
    x = x/(LastNode-1);
    node(LastNode, :) = x;
    CentNode = LastNode;
  
    for j = 1:n
        m = j+1;
        if j==n
            m=1;
        end
         tri(j,:) = [j m CentNode];
       
    end
    
if ne == n
    return;
    
else
    for l=1:100
        [noTri g] = size(tri);
        [lastNode g] = size(node);
       
        if noTri == ne
            break;
        end
        elementsFlag = 1;
        %tempNode = zeros(1,4);
        for i=1:noTri
             [cvr cvc] = size(cv);
            for j=1:3
                m = j+1;
                if j==3
                    m=1;
                end
         [cj rj] = find(cv == tri(i,j));
         if isempty(find(cv == tri(i,j)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == tri(i,m));
         if isempty(find(cv == tri(i,m)))
             cm = -1; rm = -1;
         end
 
    if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm);
        x(1) = (node(tri(i,j),1)+node(tri(i,m),1))/2;
        x(2) = curve{crow}(x(1));
         x = [x(1) x(2)];
   
   temp = find(ismember(node,x,'rows'));
   
if temp<=lastNode
   tempNode(j) = temp;
   else
    node(lastNode+1, :) = x;
    lastNode = lastNode + 1;
    cv(crow, cvc+1) = lastNode;
    tempNode(j) = lastNode;
end
     else
     x(1) = (node(tri(i,j),1)+node(tri(i,m),1))/2;
     x(2) = (node(tri(i,j),2)+node(tri(i,m),2))/2;
     x = [x(1) x(2)];
     temp = find(ismember(node,x,'rows'));

if temp<=lastNode
    tempNode(j) = temp;
   
    
else
    node(lastNode+1, :) = x;
    lastNode = lastNode + 1;
    tempNode(j) = lastNode;
end
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
%%%initial domain (before meshing)
        xmin = min(node(:,1)) ;
        ymin = min(node(:,2));
        xmax = max(node(:,1)) ;
        ymax = max(node(:,2));
        
        figure
        xlim([xmin-0.5, xmax+0.5]); ylim([ymin-0.5, ymax+0.5]);  
        hold on ; 
       
        for i=1:cvr
            cmax = max(node(cv(i,1),1), node(cv(i,2), 1));
            cmin = min(node(cv(i,1),1), node(cv(i,2), 1));
          fplot(curve{i},[cmin cmax], 'b');
          hold on;
        end
        for j = 1:n
        plot(node(j,1), node(j,2), 'o', 'MarkerSize', 4,'MarkerFaceColor', 'red');
        end
         title ('The domain to be meshed')
         xlabel('x axis');
         ylabel('y axis');
end



function [node quad cv] = mesh_BEC_2L(quad,node,  curve, cv)
[ne tem] = size(quad);
[l temp2] = size(node);
[cvr cvc] = size(cv);
for e = 1:ne
        [cj rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (node(quad(e,1),1)+node(quad(e,2),1))/2;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + node(quad(e,2),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
    
      
         [cj rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (node(quad(e,2),1)+node(quad(e,3),1))/2;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + node(quad(e,3),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
      
      [cj rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (node(quad(e,3),1)+node(quad(e,4),1))/2;
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
     temp = (node(quad(e,3), :) + node(quad(e,4),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
     end
      
        [cj rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (node(quad(e,4),1)+node(quad(e,1),1))/2;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
     
     temp = (node(quad(e,4), :) + node(quad(e,1),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
    
    node(l+1,:) = (node(quad(e,1), :) + node(quad(e,2),:) + node(quad(e,3),:)+ node(quad(e,4),:))/4;
      l = l+1;
      q5 = l;
    
     
    
    newquad(e,:) = [quad(e,1) q1 quad(e,2) q2 quad(e,3) q3 quad(e,4) q4 q5];
end
quad = newquad;
return;
end


function [node quad cv] = mesh_BEC_2S(quad,node,  curve, cv)
[ne tem] = size(quad);
[l temp2] = size(node);
[cvr cvc] = size(cv)
for e = 1:ne
        [cj rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (node(quad(e,1),1)+node(quad(e,2),1))/2;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + node(quad(e,2),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
    
      
         [cj rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (node(quad(e,2),1)+node(quad(e,3),1))/2;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + node(quad(e,3),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
      
      [cj rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (node(quad(e,3),1)+node(quad(e,4),1))/2;
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
     temp = (node(quad(e,3), :) + node(quad(e,4),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q <= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
     end
      
        [cj rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (node(quad(e,4),1)+node(quad(e,1),1))/2;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
     
     temp = (node(quad(e,4), :) + node(quad(e,1),:))/2;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
      end
    
    node(l+1,:) = (node(quad(e,1), :) + node(quad(e,2),:) + node(quad(e,3),:)+ node(quad(e,4),:))/4;
      l = l+1;
      q5 = l;
     
    
    newquad(e,:) = [quad(e,1) q1 quad(e,2) q2 quad(e,3) q3 quad(e,4) q4];
end
quad = newquad;
return;

end


function [node quad cv] = mesh_BEC_3L(quad,node,  curve, cv)

[ne tem] = size(quad);
[l temp2] = size(node);
[cvr cvc] = size(cv);
for e = 1:ne
    
        [cj, rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (2*node(quad(e,1), 1) + node(quad(e,2),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
         else
    temp = (2*node(quad(e,1), :) + node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
     end
    
       [cj, rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
           crow = intersect(cj,cm)
           x(1) = (node(quad(e,1), 1) + 2*node(quad(e,2),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + 2*node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end
    
          [cj, rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
    if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (2*node(quad(e,2), 1) + node(quad(e,3),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
     else
    temp = (2*node(quad(e,2), :) + node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
            end
            
       [cj, rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
          x(1) = (node(quad(e,2), 1) + 2*node(quad(e,3),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + 2*node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
            end      
    [cj, rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (2*node(quad(e,3), 1) + node(quad(e,4),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
        else
    temp = (2*node(quad(e,3), :) + node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
        end
      [cj, rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
          x(1) = (node(quad(e,3), 1) + 2*node(quad(e,4),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,3), :) + 2*node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end 
    
      [cj, rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (2*node(quad(e,4), 1) + node(quad(e,1),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
        else
    temp = (2*node(quad(e,4), :) + node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
        end
  [cj, rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if cj == cm
          x(1) = (node(quad(e,4), 1) + 2*node(quad(e,1),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,4), :) + 2*node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
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


function [node quad cv] = mesh_BEC_3S(quad,node,  curve, cv)
[ne tem] = size(quad);
[l temp2] = size(node);
[cvr cvc] = size(cv);
for e = 1:ne
    
                   [cj, rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (2*node(quad(e,1), 1) + node(quad(e,2),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
         else
    temp = (2*node(quad(e,1), :) + node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q1 = q;
    else
        q1 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
     end
    
       [cj, rj] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
           crow = intersect(cj,cm)
           x(1) = (node(quad(e,1), 1) + 2*node(quad(e,2),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,1), :) + 2*node(quad(e,2),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q2 = q;
    else
        q2 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end
    
          [cj, rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
    if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (2*node(quad(e,2), 1) + node(quad(e,3),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
     else
    temp = (2*node(quad(e,2), :) + node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q3 = q;
    else
        q3 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
            end
            
       [cj, rj] = find(cv == quad(e,2));
         if isempty(find(cv == quad(e,2)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
          x(1) = (node(quad(e,2), 1) + 2*node(quad(e,3),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,2), :) + 2*node(quad(e,3),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q4 = q;
    else
        q4 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
            end      
    [cj, rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (2*node(quad(e,3), 1) + node(quad(e,4),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
        else
    temp = (2*node(quad(e,3), :) + node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q5 = q;
    else
        q5 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
        end
      [cj, rj] = find(cv == quad(e,3));
         if isempty(find(cv == quad(e,3)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
          x(1) = (node(quad(e,3), 1) + 2*node(quad(e,4),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,3), :) + 2*node(quad(e,4),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q6 = q;
    else
        q6 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end 
    
      [cj, rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if isempty(intersect(cj,cm)) == 0
        crow = intersect(cj,cm)
           x(1) = (2*node(quad(e,4), 1) + node(quad(e,1),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
        else
    temp = (2*node(quad(e,4), :) + node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q7 = q;
    else
        q7 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
        end
  [cj, rj] = find(cv == quad(e,4));
         if isempty(find(cv == quad(e,4)))
             cj = 0; rj = 0;
         end
         [cm, rm] = find(cv == quad(e,1));
         if isempty(find(cv == quad(e,1)))
             cm = -1; rm = -1;
         end
     if cj == cm
          x(1) = (node(quad(e,4), 1) + 2*node(quad(e,1),1))/3;
 
           x(2) = curve{crow}(x(1));
           x = [x(1) x(2)];
           temp = x;
           q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        cv(crow,cvc+1) = l;
        node(l,:) = temp;
    end
      else
    temp = (node(quad(e,4), :) + 2*node(quad(e,1),:))/3;
    q = find(ismember(node, temp, 'rows'));
    if q<= l
        q8 = q;
    else
        q8 =  l+1;
        l = l+1;
        node(l,:) = temp;
    end
         end 
    
    newquad(e,:) = [quad(e,1) q1 q2 quad(e,2) q3 q4 quad(e,3) q5 q6 quad(e,4) q7 q8];
    
end
quad = newquad;
return;
end


function print(quad, node)
[m, n] = size(quad);
[k, l] = size(node);
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

%%%  Plot the Meshed Domain
 function drawDomain_linear(q, c)
        [k, l] = size(q);
        [r, s] = size(c);
        
 
              
        xmin = min(c(:,1)) ;
        ymin = min(c(:,2));
        xmax = max(c(:,1)) ;
        ymax = max(c(:,2));
        
        figure
        xlim([xmin-0.5, xmax+0.5]); ylim([ymin-0.5, ymax+0.5]);  
        hold on ; 
        
        for i=1:k
            cx = sum(c(q(i,:), 1))/l;
            cy = sum(c(q(i,:), 2))/l;
            if (l<16)
            for j = 1:l
                m = j+1;
                if l == 9
                x1 = c(q(i,m),1);
                y1 = c(q(i,m),2);
                
                if (abs([cx cy] - [x1 y1])<=10^-10)
                    m = 1;
                x = [c(q(i,j),1) c(q(i,m),1)];
                y = [c(q(i,j),2) c(q(i,m),2)];
                plot(x,y);
                hold on; 
                break; 
                end
                end
           
                if j==l && m ~ 1;
                    m=1;
                end
                
                x = [c(q(i,j),1) c(q(i,m),1)];
                y = [c(q(i,j),2) c(q(i,m),2)];
               
                plot(x,y);
                hold on ;
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
                hold on;
            end
             end   
       
           
           text(cx-0.06,cy,'(');
             hold on
            text(cx-0.04,cy,int2str(i));
            hold on
            text(cx-0.02,cy,')');
            hold on
            
        end
        
       for i = 1:k
            for  j = 1:l
                if (l == 4) || (l == 8) || (l == 12)
                           plot(c(q(i,j),1),c(q(i,j),2), 'o', 'MarkerSize', 4,'MarkerFaceColor', 'red');
                           hold on;
                           text(c(q(i,j),1)-0.05, c(q(i,j),2)-0.05, int2str(q(i,j)));
                           hold on;
                elseif (l == 9)
                    if (j == 9)
                        plot(c(q(i,j),1),c(q(i,j),2), 'o', 'MarkerSize', 4,'MarkerFaceColor', 'blue');
                         hold on;
                         text(c(q(i,j),1)-0.05, c(q(i,j),2)-0.05, int2str(q(i,j)));
                          hold on;
                    else
                        plot(c(q(i,j),1),c(q(i,j),2), 'o', 'MarkerSize', 4,'MarkerFaceColor', 'red');
                        hold on;
                       text(c(q(i,j),1)-0.05, c(q(i,j),2)-0.05, int2str(q(i,j)));
                        hold on;
                    end
                
                elseif (l == 16)
                        if (j>=13)
                            plot(c(q(i,j),1),c(q(i,j),2), 'o', 'MarkerSize', 4,'MarkerFaceColor', 'blue');
                             hold on;
                              text(c(q(i,j),1)-0.05, c(q(i,j),2)-0.05, int2str(q(i,j)));
                               hold on;
                        else
                             plot(c(q(i,j),1),c(q(i,j),2), 'o', 'MarkerSize', 4,'MarkerFaceColor', 'red');
                              hold on;
                              text(c(q(i,j),1)-0.05, c(q(i,j),2)-0.05, int2str(q(i,j)));
                              hold on;
                        end  
         
            
                end
            end
       end
        title (['Mesh with  ' num2str(k) '  ' num2str(l) '-noded quadrilateral elements and  ' num2str(r) '  nodes']);
        xlabel('x axis');
        ylabel('y axis');

        
     end
  
 function drawDomain_os(t, c, icv, gg)           
        [k, l] = size(t);
        [r, s] = size(c);
        [cvr, cvc] = size(icv)
        nodesInElm = l;
                 p = rem(l,3); 
                if rem(l,3) ~= 0
                    l = l - p; 
                end
     
        xmin = min(c(:,1)) ;
        ymin = min(c(:,2));
        xmax = max(c(:,1)) ;
        ymax = max(c(:,2));
        
        figure
        xlim([xmin-0.5, xmax+0.5]); ylim([ymin-0.5, ymax+0.5]);  
        hold on ; 
    
            cmax = max(c(icv(1),1), c(icv(cvc), 1));
            cmin = min(c(icv(1),1), c(icv(cvc), 1));
            fplot(gg,[cmin cmax], 'b');
      
         hold on;
        
    for i=1:k
            cx = sum(c(t(i,:), 1))/nodesInElm;
            cy = sum(c(t(i,:), 2))/nodesInElm;
            for j = 1:l
                m = j+1;
                if j == l 
                    m = 1;
                end
                x = [c(t(i,j),1) c(t(i,m),1)];
                y = [c(t(i,j),2) c(t(i,m),2)];
                  if y == gg(x)
                     continue;  
                  else
                    plot(x,y);
                    hold on;
                  end 
            end     
            if p~= 0
                  plot(c(t(i,j+p),1), c(t(i,j+p),2), 'Marker', 'O', 'MarkerFaceColor', 'red')
            end
             %text(cx-0.06,cy,'(');
            % hold on
           % text(cx+0.02,cy,int2str(i));
           % hold on
            %text(cx-0.02,cy,')');
            %hold on
            
            
        end
       
        title (['Mesh with  ' num2str(k) '  ' num2str(l) '-noded quadrilateral elements and  ' num2str(r) '  nodes'])
        xlabel('x axis');
        ylabel('y axis');
 end
  
     function drawDomain(t, c, cv, gg)
        [k, l] = size(t);
        [r, s] = size(c);
        [cvr, cvc] = size(cv)
        
 
              
        xmin = min(c(:,1)) ;
        ymin = min(c(:,2));
        xmax = max(c(:,1)) ;
        ymax = max(c(:,2));
        
        figure
        xlim([xmin-0.5, xmax+0.5]); ylim([ymin-0.5, ymax+0.5]);  
        hold on ; 
        
        for i=1:cvr
            cmax = max(c(cv(i,1),1), c(cv(i,2), 1));
            cmin = min(c(cv(i,1),1), c(cv(i,2), 1));
          fplot(gg,[cmin cmax], 'b');
        end
         hold on;
        
        for i=1:k
            cx = sum(c(t(i,:), 1))/l;
            cy = sum(c(t(i,:), 2))/l;
            
                for j = 1:l
                m = j+1;
                if j == 3 || j == 6
                m = 1;
                end
                x = [c(t(i,j),1) c(t(i,m),1)];
                y = [c(t(i,j),2) c(t(i,m),2)];
                if y==gg(x)
                 continue;  
                else
                    plot(x,y);
                end
                end
            
                
            if l == 8 || l == 9
            for j = 1:8
                m = j+1;
                if j == 8
                m = 1;
                end
                x = [c(t(i,j),1) c(t(i,m),1)];
                y = [c(t(i,j),2) c(t(i,m),2)];
                if y==gg(x)
                 continue;  
                else
                    plot(x,y);
                end
               
            end
            if (l == 9)
                j = 9;
                 plot(c(t(i,j),1),c(t(i,j),2), 'o', 'MarkerSize', 4,'MarkerFaceColor', 'red');
                           hold on
            end
            end 
           %text(cx-0.06,cy,'(');
            % hold on
            %text(cx-0.04,cy,int2str(i));
            %hold on
            %text(cx-0.02,cy,')');
            %hold on
  end
        
       
        title (['Mesh with  ' num2str(k) '  ' num2str(l) '-noded quadrilateral elements and  ' num2str(r) '  nodes'])
        xlabel('x axis');
        ylabel('y axis');

        
     end
     
      function drawDomain_DCF(t, c, cv, curve)
     
        [k l] = size(t)
        [rc sc] = size(c);
        [cvr cvc] = size(cv)
        
         nodesInElm = l;
                 p = rem(l,3); 
                if rem(l,3) ~= 0
                    l = l - p; 
                end
        
 
              
        xmin = min(c(:,1)) ;
        ymin = min(c(:,2));
        xmax = max(c(:,1)) ;
        ymax = max(c(:,2));
        
        figure
        xlim([xmin-0.5, xmax+0.5]); ylim([ymin-0.5, ymax+0.5]);  
        hold on ; 
        
        for i=1:cvr
            cmax = max(c(cv(i,1),1), c(cv(i,2), 1));
            cmin = min(c(cv(i,1),1), c(cv(i,2), 1));
          fplot(curve{i},[cmin cmax], 'b');
        end
         hold on;
        
       for i=1:k
            cx = sum(c(t(i,:), 1))/nodesInElm;
            cy = sum(c(t(i,:), 2))/nodesInElm;
           
                for j = 1:l
                m = j+1;
                if j == l
                m = 1;
                end
                x = [c(t(i,j),1) c(t(i,m),1)];
                y = [c(t(i,j),2) c(t(i,m),2)];
                for r = 1:cvr
                    if y == curve{r}(x) 
                        test(r) = 1;
                    else 
                      test(r) = 0;
                    end
                end
                
                 if isempty(find(test == 1))
                     plot(x,y);
                 else
                    % disp('no plot');
                 end
                end
                
                if p~= 0
                  plot(c(t(i,j+p),1), c(t(i,j+p),2), 'Marker', 'O', 'MarkerFaceColor', 'red')
                end
            
            %text(cx,cy,int2str(i));
            %hold on
      end
 
        title (['Mesh with  ' num2str(k) '  ' num2str(l) '-noded quadrilateral elements and  ' num2str(rc) '  nodes'])
        xlabel('x axis');
        ylabel('y axis');

        
     end
