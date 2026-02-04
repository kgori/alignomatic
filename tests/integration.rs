use std::fs;
use std::path::Path;
use std::process::Command;

#[test]
fn test_bam_input_no_index() {
    let output_dir = "test_integration_output";
    let bam_input = "data/test.sorted.bam";

    // Clean up previous run if exists
    if Path::new(output_dir).exists() {
        fs::remove_dir_all(output_dir).unwrap();
    }

    // Run the command
    // We use cargo run to ensure we are testing the current code
    // Note: This might be slower, but ensures correctness without hardcoding binary paths
    let status = Command::new("cargo")
        .args(&[
            "run",
            "--",
            "--bam-input",
            bam_input,
            "--output-folder",
            output_dir,
        ])
        .status()
        .expect("Failed to execute process");

    assert!(status.success(), "Program failed to execute successfully");

    // Check if output directory was created
    assert!(
        Path::new(output_dir).exists(),
        "Output directory was not created"
    );

    // Check for config.json presence as a basic success marker
    assert!(
        Path::new(output_dir).join("config.json").exists(),
        "config.json not found in output"
    );

    // Cleanup
    if Path::new(output_dir).exists() {
        fs::remove_dir_all(output_dir).unwrap();
    }
}
