#
#  CreateRunDirectory.jl
#  CollagenPolymer
#
#  Created by Christopher Revell on 02/05/2020.
#
#

module CreateRunDirectory

# Import Julia packages
using Dates

# Import local program modules


function createRunDirectory(nMonomers,L,nDomains,μ,kT,ϵLJ,σ,k,rₑ,Ebend,D,tMax,outputInterval,nGrid,boxSize,intrctnThrshld)

    # Create directory for run data labelled with current time.
    foldername = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")
    mkdir("output/$(foldername)")
    # Store system parameters.
    open("output/$(foldername)/conditions.txt","w") do conditionsfile
        println(conditionsfile,"nMonomers      $nMonomers      ")
        println(conditionsfile,"L              $L              ")
        println(conditionsfile,"σ              $σ              ")
        println(conditionsfile,"nDomains       $nDomains       ")
        println(conditionsfile,"nGrid          $nGrid          ")
        println(conditionsfile,"intrctnThrshld $intrctnThrshld ")
        println(conditionsfile,"μ              $μ              ")
        println(conditionsfile,"kT             $kT             ")
        println(conditionsfile,"ϵLJ            $ϵLJ            ")
        println(conditionsfile,"k              $k              ")
        println(conditionsfile,"rₑ             $rₑ             ")
        println(conditionsfile,"Ebend          $Ebend          ")
        println(conditionsfile,"D              $D              ")
        println(conditionsfile,"tMax           $tMax           ")
        println(conditionsfile,"outputInterval $outputInterval ")
        println(conditionsfile,"boxSize        $boxSize        ")
    end

    return foldername
end

export createRunDirectory

end
