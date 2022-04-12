[<img src="https://assets.signaloid.io/add-to-signaloid-cloud-logo-dark-v6.png#gh-dark-mode-only" alt="[Add to signaloid.io]" height="30">](https://signaloid.io/repositories?connect=https://github.com/signaloid/Signaloid-Demo-Metallurgy-BrownHamModel#gh-dark-mode-only)
[<img src="https://assets.signaloid.io/add-to-signaloid-cloud-logo-light-v6.png#gh-light-mode-only" alt="[Add to signaloid.io]" height="30">](https://signaloid.io/repositories?connect=https://github.com/signaloid/Signaloid-Demo-Metallurgy-BrownHamModel#gh-light-mode-only)

https://user-images.githubusercontent.com/86417/158057589-18e8915d-e5b8-40b8-87ab-fd93005ee60d.mov

<br/>
<br/>

# Precipitate Dislocation Model from Brown and Ham
This example shows how uncertainties in empirical model parameters affect the uncertainty distribution of the model's output, for a model of a physical process. The example implements the equation for a materials precipitate "cutting" dislocation model from Brown and Ham[^0] and shows how metallurgists can gain insight into the uncertainty of a model of a metal alloy's strength. The example highlights how the Signaloid C0 processor allows you to take unmodified programs and track uncertainty through them[^1], getting all the benefits that you would usually only be able to obtain from a hand-crafted (and time-consuming) Monte Carlo evaluation.

## Inputs
The inputs and their ranges are:
-	`gamma`:	0.15 to 0.25
-	`phi`:		0.30 to 0.45
-	`Rs`:		1 10^-8 to 3 10^-8
-	`G`:		6 10^10 to 8 10^10
-	`b`:		2.54 10^-9 to 2.54 10^-9 (i.e., constant)
-	`M`:		2.9 to 3.2.

The parameter `gamma` is the APB energy with units J/m^2, `phi` is the precipitate volume fraction, `Rs` is mean particle radius on plane with units m, `G` is the shear modulus with units Pa, `b` is the magnitude of the Burgers vector with units m, and `M` is the Taylor factor.

## Outputs
The output is the cutting stress, `σc` where
```
                    ⎛    _________________    ⎞
       ⎛ M ⋅ γ  ⎞   ⎜   ╱8.0 ⋅ γ ⋅ φ ⋅ Rs     ⎟
  σ  = ⎜─────── ⎟ ⋅ ⎜  ╱ ───────────────── - φ⎟
   c   ⎝2.0 ⋅ b ⎠   ⎝╲╱  π ⋅ G ⋅ pow(b, 2)    ⎠
```

## Repository Tree Structure
The repository contains three different variants of a simple program implementing the Brown and Ham model. The first variant (`v1`) has all the model parameters as point-valued numbers, with the Taylor factor (`M`) computed as the mean value of a number of empirical values. The second variant (`v2`) has the Taylor factor (`M`) as the distribution constructed directly from the empirical values, and the third variant (`v3`) has all parameters except the Burgers vector (`b`) as distributions.

```
.
├── README.md
├── v1
│   ├── README.md
│   └── src
│       ├── README.md
│       └── Brown-and-Ham-no-distributions.c
├── v2
│   ├── src
│   │   ├── Brown-and-Ham-with-only-Taylor-Factor-as-distribution.c
│   │   └── README.md
│   └── README.md
└── v3
    ├── README.md
    └── src
        ├── Brown-and-Ham-with-all-parameters-as-distributions.c
        └── README.md
```


<br/>
<br/>
<br/>

[^0]: Brown, L. M., and R. K. Ham. "Dislocation-particle interactions." Strengthening methods in crystals (1971): 9–135.
[^1]: Running this example on the Signaloid C0-Cloud processor uses less than 1% of the free monthly credits on the Signaloid C0-Cloud Free Tier plan.
