const setDataset = (source, dims, annoName, annotations, fromIndex = 0) => {
  let _datasets = []
  _datasets.push({
    dimensions: dims,
    source: source,
  })
  if (annotations.length === 0) {
    // if no annotations, all the data are annotated by a label.
    _datasets.push({
      // 这个 dataset 的 index 是 `1`。
      transform: {
        type: 'sort',
        config: { dimension: annoName, order: 'desc' },
      },
      fromDatasetIndex: fromIndex,
    })
  } else {
    for (let anno of annotations) {
      _datasets.push({
        // 这个 dataset 的 index 是 `1`。
        transform: {
          type: 'filter',
          config: { dimension: annoName, value: anno },
        },
        fromDatasetIndex: fromIndex,
      })
    }
  }
  return _datasets
}

const setBatchDataset = (source, dims, annoName, annotations, batch, batchName, fromIndex = 0, append = true) => {
  let _datasets = []
  _datasets.push({
    dimensions: dims,
    source: source,
  })
  if (batchName !== null) {
    if (annotations.length === 0) {
      // if no annotations, all the data are annotated by a label.
      _datasets.push({
        // 这个 dataset 的 index 是 `1`。
        transform: [
          {
            type: 'filter', //batchName is considered a column name in.obs
            config: { dimension: batchName, value: batch },
          },
          {
            type: 'sort',
            config: { dimension: annoName, order: 'desc' },
          },],
        fromDatasetIndex: fromIndex,
      })
    } else {
      // add datasets with batch filter
      for (let anno of annotations) {
        _datasets.push({
          transform: [
            {
              type: 'filter', //batchName is considered a column name in.obs
              config: { dimension: batchName, value: batch },
            },
            {
              type: 'filter',
              config: { dimension: annoName, value: anno },
            },
          ],
          fromDatasetIndex: fromIndex,
        })
      }
      // add datasets without batch filter for 3D view
      if (append) {
        for (let anno of annotations) {
          _datasets.push({
            transform: [
              {
                type: 'filter',
                config: { dimension: annoName, value: anno },
              }
            ],
            fromDatasetIndex: fromIndex,
          })
        }
      }
    }
  } else {
    // add datasets without batch filter
    for (let anno of annotations) {
      _datasets.push({
        transform: [
          {
            type: 'filter',
            config: { dimension: annoName, value: anno },
          }
        ],
        fromDatasetIndex: fromIndex,
      })
    }
    if (append) {  // append annotations for 3D view
      for (let anno of annotations) {
        _datasets.push({
          transform: [
            {
              type: 'filter',
              config: { dimension: annoName, value: anno },
            }
          ],
          fromDatasetIndex: fromIndex,
        })
      }
    }
  }
  return _datasets
}

const Dataset = {
  setDataset,
  setBatchDataset,
}

export default Dataset