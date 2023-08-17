// Handle file drop
function dropHandler(event) {
    event.preventDefault();

    // Retrieve the files from the drop event
    const files = event.dataTransfer.files;

    // Pass the files to the handling function
    handleFiles(files);
}

// Handle drag over event
function dragOverHandler(event) {
    event.preventDefault();
    // You may add visual effects for the drag-over event here
}

// Handle the files either dropped or selected via the input
function handleFiles(files) {
    // You can loop through the FileList and handle individual files as needed
    for (const file of files) {
        uploadFile(file);
    }
}

// Function to upload individual file
function uploadFile(file) {
    // Creating a new FormData object
    let formData = new FormData();
    formData.append('file', file);

    // Get the current project name from the URL
    let project = window.location.pathname.split('/')[2];

    // Making a POST request to upload the file
    fetch(`/p/${project}/upload`, {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        // Handle success
        console.log('File uploaded successfully:', data);
    })
    .catch(error => {
        // Handle error
        console.error('File upload failed:', error);
    });
}

// Function to allow click on drop zone to trigger the file input
document.getElementById('drop-zone').addEventListener('click', () => {
    document.getElementById('file-input').click();
});
