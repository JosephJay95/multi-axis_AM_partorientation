function[D,X]=heat_geodesics_singlesource(boundary_p, faces, source_vertex_id)

            t=avgedge(boundary_p,faces); %time step%
            t=t.^2;

            A=massmatrix(boundary_p,faces); %mass matrix
            LC= cotmatrix(boundary_p,faces);
            LHS= A - t*LC;

            %initial condition at every vertex
            u0= zeros(size(boundary_p,1),1);
            
            %heat source define
            heatSrcIdx=source_vertex_id;

            %% initial condition at heat source vertex groups

            u0(heatSrcIdx)=1;
            
            %linear equation for heat intergration (left hand division)
            u=LHS\u0;  

            %% step 2 - gradient vector field

            G= grad(boundary_p,faces);%gradient field with a size of [ (Fx3) x V ] 
            gradient_u=G*u; %gradient_u has (Fx3)x1
            gradient_u = reshape(gradient_u,size(faces,1),size(boundary_p,2));%reshaping --- tet size x 3 coords
            gradient_u= gradient_u ./ normrow(gradient_u);

            X=-gradient_u;

            %% step 3 - solve poisson equation
            % This is where the gradient vector field is converted into a
            % scalar field.

            D= div(boundary_p,faces);
            divX= D* X(:); %D*X along doesnt work. So vectorize X by X(:)
            phi=LC\divX;

            % Computing the distance difference between each vertex and source vertex
            % i.e., assigning 0 distance to source vertex
            phi=phi-phi(heatSrcIdx);
            
            %Reassigning the distance values
            D=phi;

            %making sure we get absolute value of distance (not a must)
            [~,mi] = max(abs(D));
            D = D*sign(D(mi));
end