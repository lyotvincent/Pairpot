import { SearchOutlined } from '@ant-design/icons'
import React, { useRef, useState } from 'react'
import Highlighter from 'react-highlight-words'
import { Button, Input, Space, Table, Pagination } from 'antd'
var data = [
  {
    key: 1,
    title: 'ST0002_01.h5ad',
    species: 'Homo sapiens',
    tissues: 'Skin',
    technology: '10x Visium',
    pmid: '37312623',
  },
  {
    key: 2,
    title: 'ST0002_02.h5ad',
    species: 'Homo sapiens',
    tissues: 'Brain',
    technology: '10x Visium',
    pmid: '37081484',
  },
  {
    key: 3,
    title: 'ST0002_03.h5ad',
    species: 'Mus musculus',
    tissues: 'Bone',
    technology: 'GeoMx DSP',
    pmid: '37255850',
  },
  {
    key: 4,
    title: 'ST0002_04.h5ad',
    species: 'Homo sapiens',
    tissues: 'Spinal cord',
    technology: '10x Visium',
    pmid: '37095395',
  },
]

const SampleTable = () => {
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
      dataIndex: 'key',
      key: 'key',
      width: '5%',
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
      dataIndex: 'technology',
      key: 'technology',
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
  return (
    <Table
      columns={columns}
      dataSource={data}
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
export default SampleTable
