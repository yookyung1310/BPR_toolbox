# BPR_toolbox
It provides BPR (Blink-locked Pupillary Response) correction toolbox and other prior preprocessing codes for pupillometry.

It is not well known that eyeblink disturbs pupillometry by transient constriction and re-dilation after blinking offset. We call the phenomenon BPR (Blink-locked Pupillary Response). 

The toolbox infers BPR component in your data with the probabilistic inference model, and subtract the component from your data. It would improve the statistical power of your pupil data. Moreover, it corrects BPR on a blink-by-blink basis. Therefore, it deals with the substantial blink-by-blink variability of BPR, and various patterns of non-BPR components, which is impossible with heuristic methods. 

It is easy to apply your pupil data. It requires pupil and blinks data, and sampling rate only as inputs. See bprcorrect.m, the main function that corrects BPR, for details.

The BPR correction toolbox can be applied to any pupil data if the data is preprocessed. The rest preprocessing codes here was written for pupil data from Eyelink, but you would be able to apply the preprocessing codes as well if you change the code slightly. 

Please read Yoo et al. (2021) to identify what's going on in pupillometry when a subject blinks. Also, it will introduce a model and algorithm for correction, and validate it.

We hope researchers using pupillometry benefit from the toolbox, and we always welcome researchers who want to comment, modify, or improve the toolbox.

# How to run the tutorial code

1. Please copy the BPRtoolbox_tutorial folder to your Desktop folder, and run BPRtoolbox_tutorial.m in the BPRtoolbox_tutorial folder. (Type N to the question when you run the code first. It will apply pupil data preprocessing and BPR correction both. Type Y if it has been processed, then it will load the previously preprocessed data and run BPR correction only) 

2. It will show how BPR works with 2 example subjects' data from the fixation task (Yoo et al. 2021). It will take 10-30 minutes depending on your hardware in this case. 

3. Please read BPR toolbox tutorial keynote.pdf to appreciate the preprocessing and correction result with simple descriptions.

# Reference

Yoo, K., Ahn, J., & Lee, S. H. (2021). The confounding effects of eye blinking on pupillometry, and their remedy. PloS one, 16(12), e0261463.



Last update: 2022.01.14 Kyung Yoo
 
Copyright 2021 Kyung Yoo

Contact: yookyung1310@gmail.com
