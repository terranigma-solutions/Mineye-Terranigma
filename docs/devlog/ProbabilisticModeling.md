## Brainstorming

- [x] Make a forward uncertainty model for the orientations
- [ ] Gravity inversion
  - [x] Tidy up code
  - [x] More plots
  - [ ] Run:
      - Start with num_samples=200, warmup_steps=200, num_chains=2 (minimum viable)
        Once model ish working, scale up to num_samples=1000, warmup_steps=1000, num_chains=4
        Your current settings (20 samples, 5 warmup, 1 chain) are essentially "checking if the code runs" - not doing actual inference.
- [x] Play with VI
- [x] Do sphinx
- [x] Retry graph with the new factory
- [ ] Look into the scaling
- [ ] Mapping gravity in a fancier way
- Magnetic inversion
- Double inversion 


## Notes 

- Model 1:
    - Dip 1 and 5 (and 4) seems to carry all the information