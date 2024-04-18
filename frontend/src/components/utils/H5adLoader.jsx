// this component with a .h5ad file as input, and return a json file for visualization

import * as hdf5 from 'jsfive'
const H5adLoader = (file, event) => {
  let newTitle = file.name.replace(/\.h5ad$/, '')
  var f = new hdf5.File(event.target.result, file.name)
  let _obs = f.get('obs').keys
  let _index
  if (_obs.includes('_index')) {
    _index = f.get('obs/_index').value
  } else {
    console.log(f.get('obs').attrs)
    let _indexName = f.get('obs').attrs['_index']
    _index = f.get('obs/' + _indexName).value
  }

  let _value = { index: _index }
  let _clusters
  if (_obs.includes('__categories')) {
    // old h5ad versions have __categories in obs
    // load categories in obs
    let _cateNames = f.get('obs/__categories').keys
    for (let obs of _cateNames) {
      let _cate = f.get('obs/__categories/' + obs).value
      _value[obs] = f.get('obs/' + obs).value.map((item) => {
        return _cate[item]
      })
    }
    // load values in obs
    let _valNames = f
      .get('obs')
      .keys.filter(
        (item) =>
          !(
            _cateNames.includes(item) ||
            ['__categories', '_index'].includes(item)
          )
      )
    for (let key of _valNames) {
      _value[key] = f.get('obs/' + key).value
    }
    _clusters = _cateNames.map((item, id) => ({
      value: id,
      label: item,
      attr: 'categories',
    }))
    _clusters = [
      ..._clusters,
      ..._valNames.map((item, id) => ({
        value: id + _cateNames.length,
        label: item,
        attr: 'values',
      })),
    ]
  } else {
    // new version indicates categories as groups and values as datasets
    let _type = []
    for (let key of _obs) {
      if (Array.isArray(f.get('obs/' + key).keys)) {
        let _cate = f.get('obs/' + key + '/categories').value
        _value[key] = f.get('obs/' + key + '/codes').value.map((item) => {
          return _cate[item]
        })
        _type.push('categories')
      } else {
        _value[key] = f.get('obs/' + key).value
        _type.push('values')
      }
    }
    _clusters = _obs.map((item, id) => ({
      value: id,
      label: item,
      attr: _type[id],
    }))
  }
  //h5info['clusters'] = _clusters

  let _obsmNames = f.get('obsm').keys
  console.log(_obsmNames)
  let _len = _index.length
  for (let obsm of _obsmNames) {
    let _embd = f.get('obsm/' + obsm).value
    if (typeof _embd !== 'undefined') {
      let _dims = _embd.length / _len
      let _loadDims = _dims
      if (obsm === 'X_pca') {
        _loadDims = 2
      }
      if (obsm === 'AUCell_rankings') {
        _loadDims = 2
      }
      for (let dim = 0; dim < _loadDims; dim++) {
        _value[`${obsm}_${dim}`] = _index.map((_, id) => {
          return _embd[id * _dims + dim]
        })
      }
    }

  }

  let __data = _index.map((_, id) =>
    Object.entries(_value).reduce((obj, [key, value]) => {
      obj[key] = value[id]
      return obj
    }, {})
  )
  let h5info = {
    title: newTitle,
    clusters: _clusters,
    data: __data,
    embdOps: _obsmNames.map((item, id) => ({ value: id, label: item }))
  }
  return h5info
}

export default H5adLoader