# model parameters
infectset = model.infectset
vaxset = model.vaxset

# recovery immunity arguments
spr_variant = :omicron_ba1
targ_variant = :delta
recovday = 300
vaxstatus = :full
vaxrcvd = :Pfizer
vaxday = 350

days = 300:847

recov_vec = [1.0 - cs.recoveffect(today, recovday, targ_variant, spr_variant, infectset) for today in days];


vax_vec = [1.0 - cs.vaxeffect(today, infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday, mode=:spread) for today in days];
