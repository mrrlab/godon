// checkpoint creates CheckpointIO which provides various operations with checkpoints.
package checkpoint

import (
	"encoding/json"
	"time"

	"github.com/op/go-logging"

	bolt "go.etcd.io/bbolt"
)

// log is the global logging variable.
var log = logging.MustGetLogger("checkpoint")

// DATA is key name for all parameters
var DATA = []byte("parameters")

// CheckpointData stores checkpoint data.
type CheckpointData struct {
	Parameters map[string]float64
	Likelihood float64
	Iter int
	Final bool
}

// CheckpointSaver saves checkpoints.
type CheckpointIO struct {
	db     *bolt.DB
	bucket []byte
	last   time.Time
}

// NewCheckpointIO creates a new CheckpointIO.
func NewCheckpointIO(db *bolt.DB, bucket []byte) (s *CheckpointIO) {
	s = &CheckpointIO{
		db:     db,
		bucket: bucket,
	}
	return
}

// Save saves checkpoint to the database given all the values needed.
func (s *CheckpointIO) Save(data *CheckpointData) error {
	// Even if saving fails, we do not want to run this code too often.
	s.SetNow()
	dataB, err := json.Marshal(data)
	if err != nil {
		log.Error("Error serializing checkpoint", err)
		return err
	}
	err = s.db.Update(func(tx *bolt.Tx) error {
		bucket, err := tx.CreateBucketIfNotExists(s.bucket)
		if err != nil {
			return err
		}

		err = bucket.Put(DATA, dataB)
		return err
	})
	if err != nil {
		log.Error("Error saving checkpoint", err)
	}
	return err
}

// GetParameters returns map with parameter values from checkpoint.
func (s *CheckpointIO) GetParameters() (*CheckpointData, error) {
	var data *CheckpointData
	
	err := s.db.View(func(tx *bolt.Tx) error {
		bucket := tx.Bucket(s.bucket)
		if bucket == nil {
			return nil
		}

		b := bucket.Get(DATA)

		err := json.Unmarshal(b, &data)
		return err
	})

	if err != nil {
		return nil, err
	}

	if data == nil || len(data.Parameters) == 0 {
		return nil, nil
	}

	if data.Final {
		log.Noticef("Found finished likelihood optimization checkpoint (iter=%v, lnL=%v)", data.Iter, data.Likelihood)
	} else {
		log.Noticef("Found unfinished likelihood optimization checkpoint (iter=%v, lnL=%v)", data.Iter, data.Likelihood)
	}

	return data, nil
}

// Old returns true if last checkpoint save time too long ago.
func (s *CheckpointIO) Old(seconds float64) bool {
	if time.Since(s.last).Seconds() > seconds {
		return true
	}
	return false
}

// SetNow sets last checkpoint time to now.
func (s *CheckpointIO) SetNow() {
	s.last = time.Now()
}
