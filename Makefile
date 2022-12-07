TAG = weisburd/impute_sample_metadata:latest

all: build push update_sha

build:
	docker build -f docker/Dockerfile -t $(TAG) .

push:
	docker push $(TAG)

update_sha:
	docker pull $(TAG) 2>&1 | grep Digest | cut -c 9- > docker/sha256.txt
	cat docker/sha256.txt && [ ! -z "`cat docker/sha256.txt`" ] && sed -i.bak "s/${DOCKER_IMAGE_NAME}@sha256:[^\"]*/"${DOCKER_IMAGE_NAME}@`cat docker/sha256.txt`"/"  *.py
	rm *.bak
