FROM python:3.7-alpine
COPY . /app
WORKDIR /app
RUN pip install .
CMD ["depmap_omics_upload"]
