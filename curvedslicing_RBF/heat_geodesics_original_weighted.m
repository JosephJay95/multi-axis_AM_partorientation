function[D,X]=heat_geodesics_original_weighted(V,T,base_indices, weights)

            t=avgedge(V,T); %time step
            t=t.^2;

            A=massmatrix(V,T); %mass matrix
            LC= cotmatrix(V,T);
            LHS= A - t*LC;

            %initial condition at every vertex
            u0= zeros(size(V,1),1);
            %heat source define
            heatSrcIdx=base_indices;

            %% initial condition at heat source vertex groups

            u0(heatSrcIdx)=1;
            %linear equation for heat intergration (left hand division)
            u=LHS\u0;  

            %% step 2 - gradient vector field

            G= grad(V,T);%gradient field with a size of [ (Fx3) x V ] 
            gradient_u=G*u; %gradient_u has (Fx3)x1
            gradient_u = reshape(gradient_u,size(T,1),size(V,2));%reshaping --- tet size x 3 coords
            gradient_u=gradient_u ./ normrow(gradient_u);

            X=-gradient_u .* weights;

            %% step 3 - solve poisson equation
            D= div(V,T);
            divX= D* X(:); %D*X along doesnt work. So vectorize X by X(:)
            phi=LC\divX;


            %reduce distance of each vertex from source vertices
            %so that at every source V, geodesic distance is 0.

            for i=1:size(heatSrcIdx,1)
                phi=phi-phi(heatSrcIdx (i));
            end


            %making sure we get absolute value of distance
            D=phi;
            [~,mi] = max(abs(D));
            D = D*sign(D(mi));
end
