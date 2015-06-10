
HOST=`hostname`
DIR=`pwd`
CARD=mic0

echo "launching on $HOST-$CARD"
ssh $HOST-$CARD "$DIR/mic_launch.sh $DIR"

