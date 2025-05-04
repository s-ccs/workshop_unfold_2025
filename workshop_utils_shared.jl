""" 
	simulate_eeg(;
			n_repeats = 1,
			overlap = (0.5,0.2),
			noiselevel = 0,
			continuouseffect = false,
			multichannel = false)

custom simulation function, especially written for the Unfold-Workshop :)
"""
function simulate_eeg(;
seed = 1,
n_repeats = 1,
overlap = (0.5,0.2),
noiselevel = 0,
continuouseffect = false,
twobytwo = false,
multievent = false,
multichannel = false,
sfreq = 100)

	conddict =  conditions = Dict{Any,Any}(:stimulation => ["bike", "face"])
	if twobytwo
			conddict[:size] = ["small","large"]
	end

	if multievent
		conddict[:event] = ["stimulus","fixation"]
	end
	if continuouseffect 
		conddict[:sac_amp] = range(0, 15, length = 10)
	end
	
	design =
    SingleSubjectDesign(;
       conditions = conddict,
        event_order_function = UnfoldSim.Random.shuffle,
    ) |> x -> RepeatDesign(x, n_repeats);

	p1 = LinearModelComponent(; basis = p100(;sfreq), formula = twobytwo ? @formula(0~1+stimulation*size) : @formula(0 ~ 1), β = twobytwo ? [5,-1.5,2,1] : [5]);
	
	n1 = LinearModelComponent(;basis = n170(;sfreq), formula = multievent ? @formula(0 ~ 1 + event*stimulation) : @formula(0 ~ 1 + stimulation),β = multievent ? [5, 0, 0, 3] : [5,3] );
	
	p3 = LinearModelComponent(;    basis = p300(;sfreq), formula = continuouseffect ? multievent ? @formula(0~1+event*sac_amp + sac_amp^2*event) : @formula(0 ~ 1 + sac_amp + sac_amp^2) : @formula(0~1),
	    β = continuouseffect ? multievent ? [2,3, -0.5, 0.2, 0.5, -0.2] : [2, -0.5, 0.2] : [5],
	);

	onset = UniformOnset(; offset = overlap[1]*sfreq,width = overlap[2]*sfreq);

	noise = PinkNoise(; noiselevel);


	# Combine the components in a vector
	components = [p1, n1, p3]
	if multichannel
            multichannel_lab = [
                "Right Occipital Pole",
                "Left Postcentral Gyrus",
                "Left Superior Frontal Gyrus",
            ]
        
        hart = headmodel()#type = "hartmut")
        components =
            MultichannelComponent.(components, [hart => label for label in multichannel_lab])
	end
	# Simulate data
	d,e = simulate(UnfoldSim.MersenneTwister(seed), design, components, onset, noise);

	if multichannel
		d = d[1:20:end,:]
	end
	if multievent &&continuouseffect
		e.sac_amp = Unfold.allowmissing(e.sac_amp)
		e.sac_amp[e.event .== "stimulus"] .= missing	
	end
	return d,e
end
	

begin
# helper function to detect whether something is specified correctly	
function check_response(var,name,shouldbetype)
	!isa(var,shouldbetype) ? Markdown.MD(Markdown.Admonition("warning","Variable not yet defined",[Markdown.Paragraph(["The variable ", Markdown.Code(name), " is not yet correctly defined. Should be: ", Markdown.Code(string(shouldbetype))," but was:  ",Markdown.Code(string(typeof(var)))])])) : nothing
end
macro check_response(var, shouldbetype)
           name = string(var)

           return esc(:(check_response($var, $name, $shouldbetype)))
       end
end

# helper to quickly make a tip
 my_tip(header,text) = Markdown.MD(Markdown.Admonition("info","Info: "*header,[text]))

# keep this as last, because I think it needs to be outputted
HTML("""
<style>
	 div.input{
	 align:center;
	 }
input{
                display: inline-block;
                outline: 0;
                cursor: pointer;
                border-radius: 6px;
                border: 2px solid #9558B2;
                color: #9558B2;
                background: 0 0;
                padding: 8px;
                box-shadow: rgba(0, 0, 0, 0.07) 0px 2px 4px 0px, rgba(0, 0, 0, 0.05) 0px 1px 1.5px 0px;
                font-weight: 800;
                font-size: 16px;
                
	 }
 input:hover{
                    background-color: #9558B2;
                    color: #fff;
                }

</style>

""")
