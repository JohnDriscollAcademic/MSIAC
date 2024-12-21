#push!(LOAD_PATH, ".")
#push!(LOAD_PATH, ENV["MSIAC"]) # tell Julia where to find MSIAC
#using MSIAC   #(you can also write import MSIAC and then add module functions using MSIAC.function)
import MSIAC   
using PyPlot, PyCall

# we are going to plot a kernel of degree p = 2
# we are going to generate a triangular mesh
function time_test()
      nR = 2
      e1 = zeros(nR+1)
      eL = zeros(nR+1)
      eT = zeros(nR+1)
      periodic =true;
      st  = "sine_xpy"
      mIO = "mymesh.txt"
      fIO = "mySine.txt"
      for r = 0 : nR
            p = 1
            nX = 2^r * 8;
            nY = 2^r * 8 ;
            xBounds = [0.0,2.0] ; yBounds = [0.0,2.0];

            mf = MSIAC.create_2D_mesh(4,nX,nY,xBounds,yBounds,periodic,mIO,structured=true, pert = [0.0,0.0])
            # we are going to mimmic a L2 field of degree p
            t1 = @elapsed msh = MSIAC.load_mesh(mIO, structured=true)
            ff = MSIAC.create_analytic_field([p+1,p+1],["legendre","legendre"],msh,fIO,type=st)
            # we now load our data for post-processing
            t2 = @elapsed data = MSIAC.load_data(mf,ff, modalExp="Pk");  #Re-project data using a Pk basis
          #t2 = @elapsed data = MSIAC.load_2Dfield(fIO, modalExp="Legendre", degree = p);  #Re-project data using a Pk basis
            MSIAC.update_phys!(data,["legendre","legendre"],[6,6])
            t3 = @elapsed fL = MSIAC.filter_data(data, "line")
            t4 = @elapsed fT = MSIAC.filter_data(data, "tensor")

            p1, p2 = MSIAC.get_evaluation_points(data)
            exact  = MSIAC.get_exact_data(data,st,[p1,p2])  #we are cheating a bit...
            w1, w2 = MSIAC.get_evaluation_weights(data)
            """
            HAVE TO FIXE L2 ERROR NO N DEFIND
            e1[r+1] = MSIAC.get_eL2(MSIAC.get_field(data), exact, w1, w2) ;
            eL[r+1] = MSIAC.get_eL2(fL, exact, w1, w2) ;
            eT[r+1] = MSIAC.get_eL2(fT, exact, w1, w2) ;
            """
            e1[r+1] = maximum(abs.(MSIAC.get_field(data).- exact));
            eL[r+1] = maximum(abs.(fL .- exact));
            eT[r+1] = maximum(abs.(fT .- exact)) 
            np = length(p1)*length(p2)*msh.N
            @info " Times: loading mesh $t1 & field $t2 filtering Line $t3 Tensor $t4  (points = $np)"
      end
      @info " CONVERGENCE RATES"
      for r = 2:nR+1
            @info " ERROR $(e1[r]) $(eL[r]) $(eT[r]) ORDER $(log10(e1[r-1]/e1[r]) / log10(2)) $(log10(eL[r-1]/eL[r]) / log10(2)) $(log10(eT[r-1]/eT[r]) / log10(2))  "
      end
end

tf = @elapsed time_test()

