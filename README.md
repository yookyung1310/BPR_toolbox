# BPR_toolbox
It provides BPR (Blink-locked Pupillary Response) correction toolbox and other prior preprocessing codes for pupillometry.

It is not well known that eyeblink disturbs pupillomtery by transient constriction and re-dilation after blinking offset. We call the phenomenon BPR (Blink-locked Pupillary Response). 

The toolbox infers BPR component in your data with probabilistic inference model, and subtract the component from your data. It would improve statistical power of your pupil data. Moreover, it corrects BPR on a blink-by-blink basis. Therefore, it deals with the substantial blink-by-blink variability of BPR, and various pattern of non-BPR components, which is impossible with heuristic methods. 

It is easy to apply your pupil data. It requires only pupil and blink data, and sampling rate as inputs. 

The BPR correction toolbox can be applied any pupil data if the data is preprocessed. The rest preprocessing codes here was written for pupil data from Eyelink, but you would be able to apply the preprocessing codes as well if you change the code slightly. 

Please read Yoo et al. 2021 (in revision, I will notate it as soon as it is published) to indetify what's going on in pupillometry when a subject blinks. Also, it will introduce a model and algorithm for correction, and validate it.

We hope researchers using pupillomtery benefit from the toolbox, and we always welcome researchers who want to comment, modify, or improve the toolbox.

# How to run tutorial code

1. Please copy the BPRtoolbox_tutorial folder to your Desktop folder, and run BPRtoolbox_tutorial.m in the BPRtoolbox_tutorial folder. 

2. It will show how BPR works with 2 example subjects data from the fixation task (Yoo et al. 2021, In revision).

3. Please read BPR toolbox tutorial keynote.key or BPR toolbox tutorial keynote.pdf to appreciate the preprocessing and correction result with simple descriptions.

Last update: 2021.07.02 Kyung Yoo
 
Copyright 2021 Kyung Yoo

Contatct: yookyung1310@gmail.com
