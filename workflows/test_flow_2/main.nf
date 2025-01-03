process basicExample {
  input:
  val x
  output:
  stdout
  script:
  """
  echo process job $x
  """
}

workflow {
  def num = Channel.of(1,2,3)
  results = basicExample(num)
  results.view()
}