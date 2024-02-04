import { SearchOutlined } from '@ant-design/icons'
import React, { useEffect, useRef, useState } from 'react'
import Highlighter from 'react-highlight-words'
import { Button, Input, Space, Table, Pagination } from 'antd'

const DatasetTable = ({ src, col }) => {
  const [data, setData] = useState([])
  const [searchText, setSearchText] = useState('')
  const [searchedColumn, setSearchedColumn] = useState('')
  const searchInput = useRef(null)
  const handleSearch = (selectedKeys, confirm, dataIndex) => {
    confirm()
    setSearchText(selectedKeys[0])
    setSearchedColumn(dataIndex)
  }
  const handleReset = (clearFilters) => {
    clearFilters()
    setSearchText('')
  }
  const getColumnSearchProps = (dataIndex) => ({
    filterDropdown: ({
      setSelectedKeys,
      selectedKeys,
      confirm,
      clearFilters,
      close,
    }) => (
      <div
        style={{
          padding: 8,
        }}
        onKeyDown={(e) => e.stopPropagation()}>
        <Input
          ref={searchInput}
          placeholder={`Search ${dataIndex}`}
          value={selectedKeys[0]}
          onChange={(e) =>
            setSelectedKeys(e.target.value ? [e.target.value] : [])
          }
          onPressEnter={() => handleSearch(selectedKeys, confirm, dataIndex)}
          style={{
            marginBottom: 8,
            display: 'block',
          }}
        />
        <Space>
          <Button
            type="primary"
            onClick={() => handleSearch(selectedKeys, confirm, dataIndex)}
            icon={<SearchOutlined />}
            size="small"
            style={{
              width: 90,
            }}>
            Search
          </Button>
          <Button
            onClick={() => clearFilters && handleReset(clearFilters)}
            size="small"
            style={{
              width: 90,
            }}>
            Reset
          </Button>
          <Button
            type="link"
            size="small"
            onClick={() => {
              confirm({
                closeDropdown: false,
              })
              setSearchText(selectedKeys[0])
              setSearchedColumn(dataIndex)
            }}>
            Filter
          </Button>
          <Button
            type="link"
            size="small"
            onClick={() => {
              close()
            }}>
            close
          </Button>
        </Space>
      </div>
    ),
    filterIcon: (filtered) => (
      <SearchOutlined
        style={{
          color: filtered ? '#1677ff' : undefined,
        }}
      />
    ),
    onFilter: (value, record) =>
      record[dataIndex].toString().toLowerCase().includes(value.toLowerCase()),
    onFilterDropdownOpenChange: (visible) => {
      if (visible) {
        setTimeout(() => searchInput.current?.select(), 100)
      }
    },
    render: (text) =>
      searchedColumn === dataIndex ? (
        <Highlighter
          highlightStyle={{
            backgroundColor: '#ffc069',
            padding: 0,
          }}
          searchWords={[searchText]}
          autoEscape
          textToHighlight={text ? text.toString() : ''}
        />
      ) : (
        text
      ),
  })
  const columns = [
    Table.SELECTION_COLUMN,
    {
      title: 'ID',
      dataIndex: 'id',
      key: 'id',
      width: '14%',
      sorter: (a, b) => a.key - b.key,
      sortDirections: ['descend', 'ascend'],
      fixed: 'left',
    },
    {
      title: 'Title',
      dataIndex: 'title',
      key: 'title',
      width: '40%',
      ...getColumnSearchProps('title'),
    },
    {
      title: 'Species',
      dataIndex: 'species',
      key: 'species',
      width: '20%',
      filters: [
        {
          text: 'Homo sapiens',
          value: 'Homo sapiens',
        },
        {
          text: 'Mus musculus',
          value: 'Mus musculus',
        },
      ],
      filterMode: 'tree',
      filterSearch: true,
      onFilter: (value, record) => record.species === value,
    },
    {
      title: 'Tissues',
      dataIndex: 'tissues',
      key: 'tissues',
      width: '20%',
      filters: [
        {
          text: 'Brain',
          value: 'Brain',
        },
        {
          text: 'Skin',
          value: 'Skin',
        },
        {
          text: 'Bone',
          value: 'Bone',
        },
        {
          text: 'Spinal cord',
          value: 'Spinal cord',
        },
      ],
      filterMode: 'tree',
      filterSearch: true,
      onFilter: (value, record) => record.tissues === value,
    },
    {
      title: 'Technology',
      dataIndex: 'technologies',
      key: 'technologies',
      width: '20%',
      filters: [
        {
          text: '10x Visium',
          value: '10x Visium',
        },
        {
          text: 'GeoMx DSP',
          value: 'GeoMx DSP',
        },
        {
          text: 'scRNA-seq',
          value: 'scRNA-seq',
        },
        {
          text: 'Slide-seq',
          value: 'Slide-seq',
        },
      ],
      filterMode: 'tree',
      filterSearch: true,
      onFilter: (value, record) => record.technology === value,
    },
    {
      title: 'PMID',
      dataIndex: 'pmid',
      key: 'pmid',
      ...getColumnSearchProps('pmid'),
      sorter: (a, b) => a.pmid.length - b.pmid.length,
      sortDirections: ['descend', 'ascend'],
      width: '20%',
    },
    {
      title: 'Action',
      key: 'operation',
      fixed: 'right',
      width: 120,
      render: () => <a>Visualization</a>,
    },
  ]

  useEffect(() => {
    var attrStack = {}
    for (var i = 0; i < col.length; i++) {
      var colname = col[i]
      attrStack[colname] = i
    }
    var dataStack = []
    for (var i = 0; i < src.length; i++) {
      dataStack.push({
        key: src[i][0],
        id: src[i][1],
        title: src[i][2],
        species: src[i][attrStack.species],
        tissues: src[i][attrStack.tissues],
        technologies: src[i][attrStack.technologies],
        pmid: src[i][attrStack.pmid],
      })
    }
    setData(dataStack)
  }, [src, col])

  return (
    <Table
      columns={columns}
      dataSource={data}
      size="small"
      rowSelection={{}}
      scroll={{
        x: 1500,
        y: 1000,
      }}
      pagination={{
        showSizeChanger: true,
        total: data.length,
        showQuickJumper: true,
        showTotal: (total) => `Total ${total} items`,
        position: 'bottomRight',
      }}
    />
  )
}
export default DatasetTable
