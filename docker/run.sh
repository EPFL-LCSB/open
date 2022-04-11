#!/bin/sh



docker run 	--rm -it 			\
		-v $(pwd)/work:/home/open/work \
		-v $(pwd)/..:/open	\
		open_docker



#docker run 	--rm -it 			\
#		-v $(pwd)/work:/home/open/work \
#		-v $(pwd)/..:/open	\
#		-v $(pwd)/..:/openbread	\
#		open_docker

