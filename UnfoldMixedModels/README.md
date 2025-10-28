```bash
julia -e 'using Pkg;Pkg.add(["Unfold","UnfoldMakie","MixedModels","UnfoldSim","DataFrames","DataFramesMeta","WGLMakie","PlutoUI","PlutoTeachingTools","PackageCompiler"]);using PackageCompiler; create_sysimage([:Unfold,:UnfoldSim,:UnfoldMakie,:MixedModels,:DataFrames,:DataFramesMeta,:WGLMakie,:PlutoUI,:PlutoTeachingTools],sysimage_path="/home/plutoserver/workshop_unfold_2025/jl_image")'
i
julia -J jl_image
```
```julia
using Pluto
Pluto.run(host="0.0.0.0",sysimage="jl_image",require_secret_for_open_links=false,require_secret_for_access=false,warn_about_untrusted_code=false)
```

