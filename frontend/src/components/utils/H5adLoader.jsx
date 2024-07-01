// this component with a .h5ad file as input, and return a json file for visualization

import * as hdf5 from 'jsfive'

const CpdbLoader = (f) => {

  // let dataList = dataAttr.map((item) => {
  //   let tempList = f.get(`uns/dataArray/Astrocytes/${item}`).value
  //   return tempList
  // })

  // let testttt = f.get('uns/dataArray/Astrocytes').values
  // console.log(testttt)
  // let testttt1 = f.get('uns/intArray/Astrocytes').value
  // console.log(testttt1)

  let cellType = f.get(`uns/cellType`).value
  let dataKeys = f.get(`uns/dataKeys`).value
  let intArray = {}
  let cellArray = {}
  let dataArray = {}


  for (let i = 0; i < cellType.length; i++) {
    // find temp cell
    let tempCell = cellType[i]


    // add intArray
    let oldList1 = f.get(`uns/intArray/${tempCell}`).value
    // let newList1 = oldList1.map(item => ({ "0": item }))
    // intArray[tempCell] = newList1
    intArray[tempCell] = oldList1

    // console.log(oldList1)
    // console.log(newList1)

    // add cellArray
    let oldList2 = f.get(`uns/cellArray/${tempCell}`).value
    // let newList2 = oldList2.map(item => ({ "0": item })) 
    // cellArray[tempCell] = newList2
    cellArray[tempCell] = oldList2

    // add dataArray
    // let testttt = f.get('uns/dataArray/Astrocytes').values[0].value
    // ['TENM2-FLRT3-complex-ADGRL1', 'Astrocytes-Astrocytes', '0.5095541', '1.0', '0.5095541', 
    // 'nan', '0.0', '0.0', 'nan', 'False', 'Ligand-Receptor', 'Signaling by Fibronectin']

    let dataList = dataKeys.map((item) => {
      let tempList = f.get(`uns/dataArray/${tempCell}/${item}`).value
      return tempList
    })
    dataArray[tempCell] = dataList  // 将处理后的 dataList 存储到 dataArray 中

    // dataArray without T
    // let dataAttr = f.get('uns/dataArray/Astrocytes').keys
    // dataAttr.shift()
    // // 遍历查找下属性值
    // for (let i = 0; i < cellType.length; i++) {
    //   let tempCell = cellType[i]
    //   let dataList = dataAttr.map((item) => {
    //     let tempList = f.get(`uns/dataArray/${tempCell}/${item}`).value
    //     return tempList
    //   })
    //   // console.log(dataList)
    //   // 用map转置下
    //   let tempArray = dataList[0].map((_, i) => {
    //     let obj = {}
    //     dataAttr.forEach((attr, index) => {
    //       obj[attr] = dataList[index][i]
    //     })
    //     return obj
    //   })
    //   // 添加到dataArray中
    //   dataArray[tempCell] = tempArray

    // }
  }

  return {
    'cellType': cellType,
    'dataArray': dataArray,
    'intArray': intArray,
    'cellArray': cellArray,
    'dataKeys': dataKeys,
  }

  // // step0: cellType
  // cellType = f.get('uns/intArray').keys

  // // step1: intArray,cellArray
  // for (let i = 0; i < cellType.length; i++) {
  //   let tempCell = cellType[i]

  //   // intArray
  //   let oldList = f.get(`uns/intArray/${tempCell}`).value // ["EFNB1-EPHB1","PTN-PTPRZ1"] 这种类型
  //   // "EFNB1-EPHB1" -> {"0":"EFNB1-EPHB1"}
  //   let newList = oldList.map((item) => ({ "0": item }))
  //   intArray[tempCell] = newList

  //   // cellArray
  //   let oldList1 = f.get(`uns/cellArray/${tempCell}`).value
  //   let newList1 = oldList1.map((item) => ({ "0": item }))
  //   cellArray[tempCell] = newList1
  // }

  // // dataArray
  // let dataAttr = f.get('uns/dataArray/Astrocytes').keys
  // dataAttr.shift()
  // // 遍历查找下属性值
  // for (let i = 0; i < cellType.length; i++) {
  //   let tempCell = cellType[i]
  //   let dataList = dataAttr.map((item) => {
  //     let tempList = f.get(`uns/dataArray/${tempCell}/${item}`).value
  //     return tempList
  //   })
  //   // console.log(dataList)
  //   // 用map转置下
  //   let tempArray = dataList[0].map((_, i) => {
  //     let obj = {}
  //     dataAttr.forEach((attr, index) => {
  //       obj[attr] = dataList[index][i]
  //     })
  //     return obj
  //   })
  //   // 添加到dataArray中
  //   dataArray[tempCell] = tempArray

  // }
  // console.log(dataArray)




  // 外层结构：{'Astrocytes':[],...}
  // 内层列表结构：[{"0":"EFNB1-EPHB1"}, {"0":"PTN-PTPRZ1"}..]

  // for (let i = 0; i < cellType.length; i++) {
  //   let tempCell = cellType[i]

  //   // intArray
  //   let oldList = f.get('uns/intArray/' + tempCell).value // ["EFNB1-EPHB1","PTN-PTPRZ1"] 这种类型
  //   // 把oldList中的 "EFNB1-EPHB1" -> {"0":"EFNB1-EPHB1"}
  //   let newList = []
  //   for (let j = 0; j < oldList.length; j++) {
  //     let tempDic = {}
  //     tempDic["0"] = oldList[j]
  //     newList.push(tempDic)
  //   }
  //   intArray[tempCell] = newList



  // step2:封装下dataArray
  // 外层结构也是：{'Astrocytes':[],...}
  // 内层列表结构：
  /*[
    {
      "interaction_group":"EFNB1-EPHB1",
      "celltype_group":"Astrocytes-Astrocytes",
      "scaled_means":0.9170507193,
      ...
    },
  */

  // console.log(dataAttr)
  // let rows = columns[0].map((_, rowIndex) => columns.map(column => column[rowIndex]));
}

const NetLoader = (f) => {
  // 1. categories
  let dataNet = {}
  let cateList = f.get('uns/network/categories/name').value
  let dataCate = cateList.map((item) => ({ "name": item }))
  dataNet["categories"] = dataCate

  // 2.links
  let linkKeys = f.get('uns/network/links').keys
  linkKeys.shift() //TODO: 显式地去掉_index, 避免bug，其他地方同理
  let linkList = linkKeys.map((item) => f.get(`uns/network/links/${item}`).value)

  // 用map转置下
  let dataLink = linkList[0].map((_, i) => {
    let obj = {}
    linkKeys.forEach((attr, index) => {
      obj[attr] = linkList[index][i]
    })
    return obj
  })
  dataNet["links"] = dataLink

  // 3.nodes
  let nodeKeys = f.get('uns/network/nodes').keys
  nodeKeys.shift() //TODO: 显式地去掉_index, 避免bug，其他地方同理
  let nodeList = nodeKeys.map((item) => f.get(`uns/network/nodes/${item}`).value)
  // console.log(nodeList)

  let maxValue = -Infinity
  let minValue = Infinity

  // 用map转置下
  let dataNodes = nodeList[0].map((_, i) => {
    let obj = {}
    nodeKeys.forEach((attr, index) => {
      obj[attr] = nodeList[index][i]
    })
    // return obj

    // update min/max
    if (obj.hasOwnProperty('value')) {
      let value = parseFloat(obj['value'])
      if (value > maxValue) maxValue = value
      if (value < minValue) minValue = value
    }


    if (obj.hasOwnProperty('symbolSize')) {
      obj['symbolSize'] = parseFloat(obj['symbolSize'])
    }

    return obj
  })

  let frac = maxValue - minValue
  dataNodes.forEach(node => {
    if (node.hasOwnProperty('value')) {
      if (frac == 0)
        node['symbolSize'] = 20
      else
        node['symbolSize'] = 20 * ((node['value'] - minValue) / frac + 1)
    }
  })

  dataNet["nodes"] = dataNodes
  return dataNet
}

const MarkerLoader = (f) => {
  // 1. rank_genes_groups
  let _names = f.get('uns/rank_genes_groups/names').value
  let _lfc = f.get('uns/rank_genes_groups/logfoldchanges').value
  let _pval = f.get('uns/rank_genes_groups/pvals_adj').value
  let _expr = f.get('uns/rank_genes_groups/expr').value
  let _frac = f.get('uns/rank_genes_groups/frac').value
  return {
    'names': _names,
    'lfc': _lfc,
    'pval': _pval,
    'expr': _expr,
    'frac': _frac,
  }
}

const DenLoader = (f) => {
  // dendrogram
  let _denInfo = 'uns/dendrogram/dendrogram_info/'
  let _icoord = f.get(`${_denInfo}/icoord`).value
  let _dcoord = f.get(`${_denInfo}/dcoord`).value
  let _ivl = f.get(`${_denInfo}/ivl`).value
  let _leaves = f.get(`${_denInfo}/leaves`).value
  let _cl = f.get(`${_denInfo}/color_list`).value
  let _lcl = f.get(`${_denInfo}/leaves_color_list`).value
  return {
    'icoord': _icoord,
    'dcoord': _dcoord,
    'ivl': _ivl,
    'leaves': _leaves,
    'color_list': _cl,
    'leaves_color_list': _lcl
  }
}

const H5adLoader = (file, event, loadContent = ['meta']) => {
  let newTitle = file.name.replace(/\.h5ad$/, '')
  var f = new hdf5.File(event.target.result, file.name)
  let h5info = {}

  // move the relation heatmap to CpdbLoader
  let _cpdb = {}
  let _net = {}
  let _marker = {}
  let _den = {}
  if (loadContent.includes("cpdb")) {
    _cpdb = CpdbLoader(f)
  }

  if (loadContent.includes("net")) {
    _net = NetLoader(f)
  }

  if (loadContent.includes("marker")) {
    _marker = MarkerLoader(f)
  }

  if (loadContent.includes("den")) {
    _den = DenLoader(f)
  }

  if (loadContent.includes("meta")) {
    let _obs = f.get('obs').keys
    let _index
    if (_obs.includes('_index')) {
      _index = f.get('obs/_index').value
    } else {
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
    let _len = _index.length
    for (let obsm of _obsmNames) {
      let _embd = f.get('obsm/' + obsm).value
      if (typeof _embd !== 'undefined') {
        let _dims = _embd.length / _len
        let _loadDims = _dims
        if (obsm === 'X_pca') {
          _loadDims = 2
        }
        if (obsm === 'X_pca_harmony') {
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
    h5info = {
      title: newTitle,
      clusters: _clusters,
      data: __data,
      embdOps: _obsmNames.map((item, id) => ({ value: id, label: item })),
    }
  }

  let res = {
    // add h5ad data
    ...h5info,
    // add rel-heat
    ..._cpdb,
    // add network
    ..._net,
    // add marker
    ..._marker,
    // add dendrogram
    ..._den,
  }

  return res
}

export default H5adLoader