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

// MAIN is key name for all parameters
var MAIN = []byte("main")

// CheckpointData stores checkpoint data.
type CheckpointData struct {
	Parameters map[string]float64
	Likelihood float64
	Iter       int
	Final      bool
}

// CheckpointSaver saves checkpoints.
type CheckpointIO struct {
	db   *bolt.DB
	key  []byte
	last time.Time
	seconds float64
}

// NewCheckpointIO creates a new CheckpointIO.
func NewCheckpointIO(db *bolt.DB, key []byte, seconds float64) (s *CheckpointIO) {
	s = &CheckpointIO{
		db:  db,
		key: key,
		seconds: seconds,
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
	err = SaveData(s.db, s.key, dataB)
	if err != nil {
		log.Error("Error saving checkpoint", err)
	}
	return err
}

// GetParameters returns map with parameter values from checkpoint.
func (s *CheckpointIO) GetParameters() (*CheckpointData, error) {
	var data *CheckpointData

	b, err := LoadData(s.db, s.key)

	if err != nil || b == nil {
		return nil, err
	}

	err = json.Unmarshal(b, &data)

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
func (s *CheckpointIO) Old() bool {
	if time.Since(s.last).Seconds() > s.seconds {
		return true
	}
	return false
}

// SetNow sets last checkpoint time to now.
func (s *CheckpointIO) SetNow() {
	s.last = time.Now()
}

// SaveData saves values in bolt database.
func SaveData(db *bolt.DB, key []byte, data []byte) error {
	if db == nil {
		return nil
	}
	err := db.Update(func(tx *bolt.Tx) error {
		b, err := tx.CreateBucketIfNotExists(MAIN)
		if err != nil {
			return err
		}

		err = b.Put(key, data)
		return err
	})
	return err
}

// LoadData loads data from bolt database.
func LoadData(db *bolt.DB, key []byte) ([]byte, error) {
	var data []byte
	if db == nil {
		return nil, nil
	}
	err := db.View(func(tx *bolt.Tx) error {
		b := tx.Bucket(MAIN)
		if b == nil {
			return nil
		}

		v := b.Get(key)
		if v != nil {
			data = make([]byte, 0, len(v))
			copy(data, v)
		}
		return nil
	})
	if err != nil {
		return nil, err
	}
	return data, nil
}
