# BPR_toolbox
It provides BPR (Blink-locked Pupillary Response) correction toolbox and other prior preprocessing codes for pupillometry.

It is not well known that eyeblink disturbs pupillomtery by transient constriction and re-dilation after blinking offset. We call the phenomenon BPR (Blink-locked Pupillary Response). 

The toolbox infers BPR component in your data with probabilistic inference model, and subtract the component from your data. It would improve statistical power of your pupil data. Moreover, it corrects BPR on a blink-by-blink basis. Therefore, it deals with the substantial blink-by-blink variability of BPR. 

It is easy to apply your pupil data. It requires only pupil and blink data, and sampling rate as inputs. 

The BPR correction toolbox can be applied any pupil data if the data is preprocessed. The rest preprocessing codes here was written for pupil data from Eyelink, but you would be able to apply the preprocessing codes as well if you change the code slightly. 

Please read Yoo et al. 2021 (in revision, I will notate it as soon as it is published) to indetify what's going on in pupillometry when a subject blinks. Also, it will introduce a model and algorithm for correction, and validate it.

We hope researchers using pupillomtery benefit from the toolbox, and we always welcome researchers who want to comment, modify, or improve the toolbox.

Last update: 2021.06.29 Kyung Yoo
 
Copyright 2021 Kyung Yoo

Contatct: yookyung1310@gmail.com
