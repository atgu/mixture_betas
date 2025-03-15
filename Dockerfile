FROM python:3.8-slim-bullseye

RUN apt update && apt install -y wget curl python3-pip locales unzip
RUN pip3 install pandas wget numpy scipy tqdm diptest 
RUN pip install fsspec gcsfs

LABEL added_three_component=True
COPY run_optimizer_alt_ss_w_mixed_effects.py /home/
COPY run_optimizer_SE_w_mixed_effects.py /home/
COPY assign_EF.py /home/

#COPY splicing_stats /home/

