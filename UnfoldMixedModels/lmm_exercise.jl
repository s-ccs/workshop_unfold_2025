### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 9a051f09-aa5b-41ec-bb77-7ec22978bd9d
using Pkg

# ╔═╡ 56e45727-ea6d-434a-8515-cd6e60bda3a6
Pkg.activate(expanduser("~/workshop_cuttingGarden2023/sysimage"))

# ╔═╡ d9912a4c-5d3a-11ee-381e-03ad95d59994
begin
	
	
	using Unfold # LMM analysis
	using UnfoldSim # Simulation
	using UnfoldMakie,CairoMakie # Plotting + Backend
	using MixedModels # LMM backend
	using DisplayAs # better display
	using Random,StableRNGs,DataFrames,StatsBase, Distributions # some helpers
	using PlutoUI,PlutoTables, PlutoTeachingTools # better Pluto Displays
	using MakieThemes
	set_theme!(ggthemr(:fresh))
end

# ╔═╡ f65f535c-e9f3-4972-ad20-42eeaac43427
md"""
# Task 1: Get to know Pluto.jl
Interactive get-to-know julia + LMMs + EEG notebook for the **Cutting-Garden 2023**

Author: [Benedikt Ehinger](www.s-ccs.de)

---
"""

# ╔═╡ 363cb188-97c7-4ece-b416-08256da0b5f9
md"""
What follows is some simple code to get an interesting plot. The idea is to get to know Pluto.jl a little bit & implement a simple interactive slider.
"""

# ╔═╡ 47627344-61e0-4992-84a9-612a5e6798f3
PlutoTeachingTools.warning_box(md"No need to understand the code, you can simply skip it for now and look back at it at a more relaxed point in time")

# ╔═╡ 4de80b12-59c8-4cee-86d7-d13463fa263a
function lorenz_step!(state,fixed,Δt)
    x, y, z = state   # current state
	σ, ρ, β = fixed   # current coefficients
    
    dx = σ*(y-x)
    dy = x*(ρ-z) - y
    dz = x*y - β*z
	state .= [x + Δt*dx,
			y + Δt*dy,
			z + Δt*dz]
	return  copy(state)
end;

# ╔═╡ e4b321f8-e5e4-4459-aa89-9d8d9c640d4d
aside(tip(md"""
This function calculates the Lorenz-Attractor.

Notice the `!` at `lorenz_step!`, it indicates that the function is modifying a container of the input (in this case it updates `state`)
"""),v_offset=-250)

# ╔═╡ 3cbe1c9d-d3e3-4757-ba8d-2867e37b1dbe
tip(md"""
`Pluto.jl` puts the outputs on top of the cell, not below
""")

# ╔═╡ 65befddb-2542-473e-8b3b-0b5c5b9fe839
question_box(md"""Try changing one of the values below of the `parameter` Array - does the plot update?

Tip: `Shift`+`Enter` or `Ctrl`+`S` will automatically run a cell
""")

# ╔═╡ a851ef86-c830-4de3-b485-8997f9e5b81c
	parameters = [1.0,10,7/8]

# ╔═╡ f2ea0255-720b-4689-ac7d-798ef1d0c990
begin
	state0 = [1.0,0,0]
	Δt = 0.05
	res = hcat([lorenz_step!(state0,parameters,Δt) for s in 1:1000]...)
end

# ╔═╡ ae74ba42-1aaa-4498-84e2-7633c64a3a37
let # enforces local scope, you cannot access variables from here outside this cell
	if any(isnan,res) | !all(isfinite,res)
		danger(md"Can't show the figure: Your parameters created values that are Infinite - please use different parameters!")
	else
	f,ax,l = lines(res[1,:],res[2,:]; color = res[3,:])
	xlims!(ax,-8,8)
	f
	end
end

# ╔═╡ 8d0da37e-25bf-42e0-82cc-72cd24a5c82d
md"""
We can easily use Sliders instead of fixing the parameters. 

A slider is defined like this:
```julia
@bind yourvariable PlutoUI.Slider(from:1:to,show_value=true,
					default=defaultvalue) # be sure to choose something else than 0!
 
```
"""

# ╔═╡ c9f67a9a-7862-4f51-93c9-ae644fc836dc
# add your slider here

# ╔═╡ 0b009035-d91a-4c46-a762-3ae33e5bae18
question_box(md"""
1. Generate up to three sliders for the three parameters in `parameters`.
2. Remember to change the `parameters` vector by replacing the current value with the `variablename` you gave the slider!

If you want to be fancy, you can get the σ, ρ, β characters by typing e.g. `\beta` + `TAB`
""")

# ╔═╡ bd2bfaab-86f5-4914-8825-87894b6a182f
Markdown.MD(Markdown.Admonition("tip","Bonus-Question",[md"""
If you have time, provide some `PlutoUI.CheckBox` or `PlutoUI.Select` elements, to change which dimension is plotted on the x/y axis
"""]))


# ╔═╡ ea71ee25-35bf-49bb-a869-db2cfe98c6f3
md"""
# Task 2: Simulating Multi-subject EEG
"""

# ╔═╡ 4a23c228-9494-4a59-8c3f-0b4d1621e322
md"""
### Experimental Design
We start by simulating a design with a 2-level factor `condition` (🚗 vs. 😊) and one `continuous` effect from -5:5. We simulate 20 subjects, and 40 items for now.
"""

# ╔═╡ e40bf805-8479-41a6-9083-5d1471a2bc5d
begin
design = MultiSubjectDesign(;
		n_subjects = 20, n_items = 40,
        items_between =Dict(:condition=>["car","is_face? 😊"], 			
							:continuous=>range(-5,5,length=10))) 
first(generate(design),5)
end

# ╔═╡ 8371a599-2165-477c-ace3-4f3eb77b5004
md"""
### ERP effects
Next we simulating three components, a P1, a N1, and a P300. For each we have to specify the `fixed` effects `β` and the respective `random` effects `σs`.
"""

# ╔═╡ 56054c8b-02e2-4fb0-a194-fa0c562efa9c
md"""
### Simulate the data
We add some PinkNoise with a certain noiselevel (you can ignore the UniformOnset, that is for continuous data modelling)
"""

# ╔═╡ 53ad0364-f7ad-4b27-90f8-f4e06bc26c22
md"""
### Analyze the data
Let's run a 2-stage ERP analysis, extracting the intercept (condition = 🚗) and difference of condition (😊 - 🚗).
"""

# ╔═╡ cf5c8e9d-44bd-4701-910e-232972b0195d
md"""
To get an overview, here you can see the output we get - a nice 🧹tidy dataframe.
"""

# ╔═╡ cdc4b2d3-4d9d-4b56-8153-7175cd86acc4
md"""
Finally, let's plot what we have simulated.
"""

# ╔═╡ 1f9c0f4d-feff-4b92-b5f5-e03b92ff8bba
nl_slider = @bind noiselevel PlutoUI.Slider(1:2:40,show_value=true);

# ╔═╡ 82ef9822-c2db-4590-8aa2-6d8bb38264a5
σs_n1_intercept_slider = @bind σs_n1_intercept PlutoUI.Slider(0:1:10,default=1,show_value=true);

# ╔═╡ b3122077-dfae-4117-a90f-1f837a1d1196
   begin
	   
   sfreq = 50 # please keep this low to not take too much CPU of other users :-)
   p1 = MixedModelComponent(;basis=p100(;sfreq=sfreq), 
	   formula=@formula(0~1+(1|subject)+(1|item)),
	   β=[10],
	   σs=Dict(:subject=>3,:item=>[1]),contrasts=Dict())
                   
   n1 = MixedModelComponent(;basis=n170(;sfreq=sfreq),
   # note: n170 is already negative, beta doesnt need to be neg
	   formula=@formula(0~1+condition+ (1+condition|subject)+ (1+condition|item)), 
	   β=[5,-3],
	   σs=Dict(:subject=>[σs_n1_intercept,4],:item=>[0.5,0.5]),contrasts=Dict())
	   
	p3 = MixedModelComponent(;basis=p300(;sfreq=sfreq),
		formula=@formula(0~1+continuous+(1+continuous|subject)+(1+continuous|item)),
		β=[5,1],
		σs=Dict(:subject=>[4,0],:item=>[0.5,0.5]),contrasts=Dict())
   end;

# ╔═╡ 55c04898-9783-4e74-82a9-c38f66df106e
dat,evts = UnfoldSim.simulate(MersenneTwister(1),design,[p1,n1,p3],UniformOnset(1,1),PinkNoise(noiselevel=noiselevel);return_epoched=true);

# ╔═╡ 57d0f9dd-4f52-4933-9386-ab74373b76e8
md"""
Some sliders to modify the simulation

|||
|---|---|
|Noiselevel | $nl_slider| 
|N1 intercept random effect σs|$σs_n1_intercept_slider

"""

# ╔═╡ 2bde3123-09a5-4faa-84ba-8ef579179511
question_box(md"""
Looking at the single subject ERPs: At what periods & effects do you have high between-subject variability?
""")

# ╔═╡ 1fb9316e-0943-4ba9-964d-9ac8daf44eba
md"""
# Task 3: ROI based LMMs
"""

# ╔═╡ 8b91f1ba-8667-489f-a396-6da5a238cbe8
md"""
For starters we want to fit a single LMM to a Region of Interest (ROI). Given we simulated a single channel right now, we only need to extract the activity in the right time-window.
"""

# ╔═╡ 30d27b0f-9219-49fb-9da4-48c111f7676e
md"""
Calculate ROI-value, baseline, and fit the model
"""

# ╔═╡ bac2de36-f031-4033-8fe8-87ad66250491
aside(tip(md"Careful! These p-values can be missleading, they are in most cases too small!"),v_offset=-170)

# ╔═╡ e60e3222-a563-401a-ac71-7cc86a6a5972
question_box(md"""
Do we have a "significant" face effect? How large is it compared to what we simulated? Why?
""")

# ╔═╡ d78b7df1-bd30-4bca-a5c9-40dfad22a2d9
question_box(md"""
Some **Bonus** tasks:
- Add an `item` effect (for technical reasons, the item effect needs to be inserted before the subject effect)
- Use `roi_bsl` instead of `bsl` - what parameter changes? 
""")

# ╔═╡ 5e3c5c35-6ce1-4196-bda3-fdf5426650a1
md"""
For the roi-baseline task a hint:
"""

# ╔═╡ 33eb0e21-2cac-481a-bec7-0dec8f4e5fa2
hint(md"""
A checkbox might be helpful:

`$(@bind dobaseline PlutoUI.CheckBox()`

And a julia shorthand for if/else: `dobaseline ? data_bsl : data`

""")

# ╔═╡ 668ab2c9-b3bb-4354-841c-4c7116831df1
md"""
# Task 4: Mass univariate MixedModels
"""

# ╔═╡ b6b7e722-33b0-48b7-98e5-85efe4c049f6
times = range(-0.1,0.5,length=size(dat,1)); # timing vector for plotting

# ╔═╡ 5bff1caf-ac56-4ee7-8164-43b23e523e2e
begin
erp2stage = []
for s = unique(evts.subject)
	ix = evts.subject .== s # get data of single subjects
	d = @view dat[:,ix]
	
	m = fit(UnfoldModel,@formula(0~1+condition),evts[ix,:],d,times) # fit model
	r = coeftable(m) # extract results
	r.subject .= s # add subject label
	append!(erp2stage,[r])
end
erp2stage = vcat(erp2stage...)
end;

# ╔═╡ 89ba4d19-e805-41ae-81ea-81d05a2444a2
first(erp2stage,5)

# ╔═╡ 91227813-c5c1-4117-8d56-374e8796c475
describe(erp2stage)

# ╔═╡ 25f0fdf0-8ca9-4177-a310-8d831288a972
let
f = plot_erp(erp2stage;mapping=(;group=:subject))
	
# add GrandAverage ERP lines
GAerp = combine(groupby(erp2stage,[:time,:coefname]),:estimate => mean=>:estimate)
[lines!(current_axis(),gp.time,gp.estimate,color=:black,linewidth=5) for gp in groupby(GAerp,:coefname)]

f
end

# ╔═╡ c938df0c-944f-44a2-bff8-aa89ace5c95b
begin
bslwindow = times.<0
roiwindow = times .> 0.15 .&& times .< 0.2
end

# ╔═╡ f6cdd050-2283-47d9-bd21-f61474c25d8b
begin
evts.roi .= NaN
evts.bsl .= NaN

winsorizedmean(x) =  mean(winsor(x))
	
for s = unique(evts.subject)
	ix = evts.subject .== s # get data of single subjects
	d = dat[:,ix] 
	
	evts.bsl[ix] = winsorizedmean(eachrow(d[bslwindow,:])) # bsl correction
	
	evts.roi[ix] = winsorizedmean(eachrow(d[roiwindow,:]))

end
evts.roi_bsl = evts.roi .- evts.bsl # bsl correct ROI


# Fit the model!
m_roi = fit(MixedModel,@formula(roi~1+condition+(1+condition|subject)),evts);
end;

# ╔═╡ c92fbbeb-b1f9-4b3f-ba40-3d0bc88ba3cd
DisplayAs.Text(m_roi)

# ╔═╡ 980e21c2-c7c3-46fd-8d61-f6c71e79a3a7
dat_bsl = dat .- mean(dat[times.<-0.05,:],dims=1); #baseline

# ╔═╡ 06de74ec-c56b-4439-9c0d-5f04823e1a47
md"""
Baseline Correct? $(@bind dobaseline PlutoUI.CheckBox())
"""

# ╔═╡ 8315f08c-3c76-433a-b9ee-4669a1b715d0
# ╠═╡ show_logs = false
uf =fit(UnfoldModel,@formula(0~1+condition+(1|subject)),evts,dobaseline ? dat_bsl : dat,times);

# ╔═╡ 2806dd13-a328-4d43-957b-51ec11f3d0d8
plot_erp(coeftable(uf);mapping=(;col=:group),axis=(;limits=(-0.1,0.5,-6,12)))

# ╔═╡ 1c513c17-99ae-4b0d-8caf-a3c5c4d81cf3
question_box(
	md"""
	Play around with the model:
	1. Move the sliders from the top and play around with noise / subject-variability
	1. add a random slope for `condition`
	1. add an `(1|item)` effect (be sure to put it before the subject effect, else you get an error)
	"""
)

# ╔═╡ fc4bc561-ff65-414a-80fb-d2c2e3083656
TableOfContents()

# ╔═╡ Cell order:
# ╠═9a051f09-aa5b-41ec-bb77-7ec22978bd9d
# ╠═56e45727-ea6d-434a-8515-cd6e60bda3a6
# ╠═d9912a4c-5d3a-11ee-381e-03ad95d59994
# ╟─f65f535c-e9f3-4972-ad20-42eeaac43427
# ╟─363cb188-97c7-4ece-b416-08256da0b5f9
# ╟─47627344-61e0-4992-84a9-612a5e6798f3
# ╠═4de80b12-59c8-4cee-86d7-d13463fa263a
# ╟─e4b321f8-e5e4-4459-aa89-9d8d9c640d4d
# ╟─3cbe1c9d-d3e3-4757-ba8d-2867e37b1dbe
# ╠═f2ea0255-720b-4689-ac7d-798ef1d0c990
# ╟─65befddb-2542-473e-8b3b-0b5c5b9fe839
# ╠═a851ef86-c830-4de3-b485-8997f9e5b81c
# ╠═ae74ba42-1aaa-4498-84e2-7633c64a3a37
# ╟─8d0da37e-25bf-42e0-82cc-72cd24a5c82d
# ╠═c9f67a9a-7862-4f51-93c9-ae644fc836dc
# ╟─0b009035-d91a-4c46-a762-3ae33e5bae18
# ╟─bd2bfaab-86f5-4914-8825-87894b6a182f
# ╟─ea71ee25-35bf-49bb-a869-db2cfe98c6f3
# ╟─4a23c228-9494-4a59-8c3f-0b4d1621e322
# ╠═e40bf805-8479-41a6-9083-5d1471a2bc5d
# ╟─8371a599-2165-477c-ace3-4f3eb77b5004
# ╠═b3122077-dfae-4117-a90f-1f837a1d1196
# ╟─56054c8b-02e2-4fb0-a194-fa0c562efa9c
# ╠═55c04898-9783-4e74-82a9-c38f66df106e
# ╟─53ad0364-f7ad-4b27-90f8-f4e06bc26c22
# ╠═5bff1caf-ac56-4ee7-8164-43b23e523e2e
# ╟─cf5c8e9d-44bd-4701-910e-232972b0195d
# ╠═89ba4d19-e805-41ae-81ea-81d05a2444a2
# ╠═91227813-c5c1-4117-8d56-374e8796c475
# ╟─cdc4b2d3-4d9d-4b56-8153-7175cd86acc4
# ╠═25f0fdf0-8ca9-4177-a310-8d831288a972
# ╟─57d0f9dd-4f52-4933-9386-ab74373b76e8
# ╠═1f9c0f4d-feff-4b92-b5f5-e03b92ff8bba
# ╠═82ef9822-c2db-4590-8aa2-6d8bb38264a5
# ╟─2bde3123-09a5-4faa-84ba-8ef579179511
# ╟─1fb9316e-0943-4ba9-964d-9ac8daf44eba
# ╟─8b91f1ba-8667-489f-a396-6da5a238cbe8
# ╠═c938df0c-944f-44a2-bff8-aa89ace5c95b
# ╟─30d27b0f-9219-49fb-9da4-48c111f7676e
# ╠═f6cdd050-2283-47d9-bd21-f61474c25d8b
# ╠═c92fbbeb-b1f9-4b3f-ba40-3d0bc88ba3cd
# ╟─bac2de36-f031-4033-8fe8-87ad66250491
# ╟─e60e3222-a563-401a-ac71-7cc86a6a5972
# ╟─d78b7df1-bd30-4bca-a5c9-40dfad22a2d9
# ╟─5e3c5c35-6ce1-4196-bda3-fdf5426650a1
# ╟─33eb0e21-2cac-481a-bec7-0dec8f4e5fa2
# ╟─668ab2c9-b3bb-4354-841c-4c7116831df1
# ╠═b6b7e722-33b0-48b7-98e5-85efe4c049f6
# ╠═980e21c2-c7c3-46fd-8d61-f6c71e79a3a7
# ╠═8315f08c-3c76-433a-b9ee-4669a1b715d0
# ╠═2806dd13-a328-4d43-957b-51ec11f3d0d8
# ╟─06de74ec-c56b-4439-9c0d-5f04823e1a47
# ╟─1c513c17-99ae-4b0d-8caf-a3c5c4d81cf3
# ╟─fc4bc561-ff65-414a-80fb-d2c2e3083656
