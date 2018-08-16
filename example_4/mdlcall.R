#-- Executing the modelling function. Takes - coef and other arguments
#-- specified in model.R, outputs - $sol, $sol.perf and $coef.new as an mdl object.
mdl <- MDLCore(coef=coef, 
               move.step=0.002,
               state=state,
               model=model,
               target=target) 