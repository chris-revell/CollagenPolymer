#
#  createRunDirectory.jl
#  collagen-model
#
#  Created by Christopher Revell on 02/05/2020.
#
#
__precompile__()
module CreateRunDirectory

using Dates

function createRunDirectory(Nmonomers,L,a,Ndomains,μ,kT,ϵLJ,σ,k,re,Ebend,D,tmax,dt,outputInterval)

    foldername = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")
    mkdir("output/$(foldername)")
    open("output/$(foldername)/conditions.txt","w") do conditionsfile
        println(conditionsfile,"Nmonomers      $Nmonomers      ")
        println(conditionsfile,"L              $L              ")
        println(conditionsfile,"a              $a              ")
        println(conditionsfile,"Ndomains       $Ndomains       ")
        println(conditionsfile,"μ              $μ              ")
        println(conditionsfile,"kT             $kT             ")
        println(conditionsfile,"ϵLJ            $ϵLJ            ")
        println(conditionsfile,"σ              $σ              ")
        println(conditionsfile,"k              $k              ")
        println(conditionsfile,"re             $re             ")
        println(conditionsfile,"Ebend          $Ebend          ")
        println(conditionsfile,"D              $D              ")
        println(conditionsfile,"tmax           $tmax           ")
        println(conditionsfile,"dt             $dt             ")
        println(conditionsfile,"outputInterval $outputInterval ")
    end

    return foldername
end

export createRunDirectory

end
