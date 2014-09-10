package optimize

import (
	"os"
)

type Optimizer interface {
	SetOptimizable(Optimizable)
	WatchSignals(...os.Signal)
	SetReportPeriod(period int)
	Run(iterations int)
}
