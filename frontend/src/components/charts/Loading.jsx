const enterLoading = (index, setLoadings) => {
  setLoadings((prevLoadings) => {
    let newLoadings = [...prevLoadings]
    newLoadings[index] = true
    return newLoadings
  })
  setTimeout(() => {
    setLoadings((prevLoadings) => {
      let newLoadings = [...prevLoadings]
      newLoadings[index] = false
      return newLoadings
    })
  }, 60000)
}

const quitLoading = (index, setLoadings) => {
  setLoadings((prevLoadings) => {
    let newLoadings = [...prevLoadings]
    newLoadings[index] = false
    return newLoadings
  })
}

const Loading = {
  enterLoading,
  quitLoading,
}

export default Loading