.PHONY: build bash

build:
	docker build -t blastn_py:latest .

bash:
	docker run -v "./scripts/:/app/scripts/" -v "./data/:/app/data/" -it --rm blastn_py:latest

remove:
	docker rmi -f blastn_py:latest

prune:
	docker system prune -a -f
